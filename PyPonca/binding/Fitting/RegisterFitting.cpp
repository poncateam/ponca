/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <algorithm>
#include <span>
#include <string>

#include "FittingList.h"
#include "../SpatialPartitioning/pykdtree.h"

// This file is the main binding code for the Fitting module.
//
// The philosophy is to avoid lambdas as much as possible and rely
// on plain functors instead. Lambdas may or may not be supported by
// some device code. In CUDA for instance, only one layer of lambda
// is supported. Function pointers can not be used either, as host
// pointers can not be called from device code.
//
// Functors only see the views of the point clouds (see view()), which
// are the same on every device.

/**
 * \brief Adds a computation to the PyComputeObject
 *
 * The main purpose of this function is to fill the outputDimension
 * which the PyCo can not do (it only sees the id);
 *
 * This function may throw if the computation is not supported by the object.
 *
 * \param self The instance of the object
 * \param id The computation to be performed
 * \param data The input array.
 */
template <typename PyCo>
void addComputation(PyCo& self, Computation id, nb::ndarray<typename PyCo::Scalar> data)
{
    ComputationDescriptor<typename PyCo::Scalar> descriptor;
    descriptor.id         = static_cast<size_t>(id);
    descriptor.outputDims = ComputeOutputDimension<PyCo>(id, nb::ndarray<>(data));
    descriptor.inputData  = data;

    self.addComputation(std::move(descriptor));
}

/**
 * \brief Perform the compute method on a point cloud
 */
struct PerformRawCloudComputation
{
    /**
     * \tparam Co The compute obejct to perform the compute method on
     *
     * \param co The compute object
     * \param cloud The view of the point cloud
     * \param i Index into the filter list
     */
    template <typename Co, typename Cloud>
    PONCA_MULTIARCH Ponca::FIT_RESULT operator()(Co& co, const Cloud& cloud, unsigned int i) const
    {
        return co.compute(cloud.begin(), cloud.end());
    }
};

/**
 * \brief Perform the compute method on a KDTree
 */
struct PerformKDTreeComputation
{
    /**
     * Neighbors are not gathered, as it would require dynamic allocations on the device:
     * the query is run again on each pass instead.
     *
     * \tparam Co The compute obejct to perform the compute method on
     *
     * \param co The compute object
     * \param tree The view of the kdtree
     * \param i Index into the filter list
     */
    template <typename Co, typename Tree>
    PONCA_MULTIARCH Ponca::FIT_RESULT operator()(Co& co, const Tree& tree, unsigned int i) const
    {
        const auto& nf = co.getNeighborFilter();
        auto neighbors = tree.rangeNeighbors(nf.frame().center(), nf.evalScale());

        if constexpr (SupportsDevice<Co>)
            return co.computeWithIds(neighbors, tree.points());
        else
        {
            // Host-only objects (such as CNC) expect STL-like containers
            std::vector<int> indices;
            std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(indices));
            return co.computeWithIds(indices, std::span(tree.points(), tree.pointCount()));
        }
    }
};

/**
 * \brief Functor calling ExtractComputation
 *
 * Unlike a pointer to ExtractComputation, it can be called from device code.
 */
struct ExtractComputationFunctor
{
    template <typename CO, typename Scalar>
    PONCA_MULTIARCH void operator()(size_t id, CO& object, size_t i, const DeviceArrayView<Scalar>& in,
                                    DeviceArrayView<Scalar>& out) const
    {
        ExtractComputation<CO>(id, object, i, in, out);
    }
};

/**
 * \brief Performs the computation
 *
 * \param self The instance of the object
 * \param cloud The point cloud
 */
template <typename PyCo, typename Cloud>
auto computeRawPointCloud(PyCo& self, const Cloud& cloud)
{
    return self.compute(cloud, PerformRawCloudComputation{}, ExtractComputationFunctor{});
}

/**
 * \brief Performs the computation
 *
 * \param self The instance of the object
 * \param cloud The point cloud
 */
template <typename PyCo, typename PyKDTree>
auto computeKDTree(PyCo& self, const PyKDTree& cloud)
{
    return self.compute(cloud, PerformKDTreeComputation{}, ExtractComputationFunctor{});
}

/**
 * \brief Binds all compute object given by factories
 *
 * \tparam P The point type. Should be compatible with PyPoncaPointCloud
 * \tparam NF The neighbor filter
 * \tparam Diff The type of differentiation
 */
template <typename PointCloud, typename _NF, Ponca::DiffType Diff>
void RegisterComputeObjects(nb::module_& m, std::set<std::string>& list)
{
    using P       = typename PointCloud::Point;
    using Scalar  = typename P::Scalar;
    using NF      = typename _NF::NF;
    using Factory = Ponca::Factory<P, NF, Diff>;

    // List of types that are not directly templated on filters
    const std::vector<std::string> nonTemplatedFilterTypes = {"CNC"};

#ifdef __CUDACC__
    const bool cudaAvailable = true;
#else
    const bool cudaAvailable = false;
#endif

    // General properties
    Factory::foreach ([&](const auto& x) {
        using T    = decltype(x.object);
        using PyCo = PyComputeObject<T, SupportsDevice<T>>;

        std::string basename = x.name;
        basename.erase(std::remove(basename.begin(), basename.end(), ' '), basename.end());

        // Search if the current type is "not templated on filter".
        bool found = false;
        for (const auto& mark : nonTemplatedFilterTypes)
            found = found || (basename.find(mark) != std::string::npos);

        // Check for existing types. It turns out some of the type within the factory
        // may not be templated on _NF, which leads to duplicates.
        // For now, this is handled here:
        nb::handle existing = nb::type<PyCo>();
        if (existing.is_valid())
        {
            // Let's hope that code whose code goes into the if-statement never reaches
            // master branch and the user :)
            if (!found)
                throw std::runtime_error("Type for " + basename + " algready registerd (" +
                                         std::string(typeid(T).name()) + ")");
            else
                return; // Return early to avoid registering the new type again
        }

        const std::string newname = basename + PointCloud::PointName + (found ? "" : _NF::name);
        auto pyco                 = nb::class_<PyCo>(m, newname.c_str());
        pyco.def(nb::init<>());
        pyco.def("setNeighborFilter", &PyCo::setNeighborFilter, nb::arg("centers"), nb::arg("radius"),
                 nb::arg("normal").none());
        pyco.def("addComputation", &addComputation<PyCo>, nb::arg("id"), nb::arg("data").none());
        pyco.def("compute", &computeRawPointCloud<PyCo, PointCloud>);
        pyco.def("compute", &computeKDTree<PyCo, PyKDTree<PointCloud>>);
        pyco.attr("isCudaAvailable") = SupportsDevice<T> && cudaAvailable;

        list.insert(basename);
    });
}

template <typename Scalar, unsigned int Dim, template <class> class NF>
void RegisterComputeObjects(nb::module_& m, std::set<std::string>& list)
{
    using namespace Ponca;
    using PointCloud = PyPointCloud<Scalar, Dim>;
    using Point      = typename PointCloud::Point;

    RegisterComputeObjects<PointCloud, NF<Point>, Ponca::FitSpaceDer>(m, list);
}

template <typename Scalar, unsigned int Dim>
void RegisterComputeObjects(nb::module_& m, std::set<std::string>& list)
{
    RegisterComputeObjects<Scalar, Dim, SWFilter>(m, list);
    RegisterComputeObjects<Scalar, Dim, CWFilter>(m, list);
    RegisterComputeObjects<Scalar, Dim, NWFilter>(m, list);
}

template <typename Scalar>
void RegisterComputationResult(nb::module_& m)
{
    const std::string computationResultName = "ComputationResults" + MangleType<Scalar>();

    auto result = nb::class_<ComputationResults<Scalar>>(m, computationResultName.c_str());
    result.def_rw("fitStatus", &ComputationResults<Scalar>::fitStatus, nb::rv_policy::reference);
    result.def_rw("results", &ComputationResults<Scalar>::results, nb::rv_policy::reference);
}

/**
 * \brief Register all instances of ComputeObjects
 *
 * \param m The module to register instances within
 */
void RegisterFitting(nb::module_& m, nb::module_& internal)
{
    // Due to registry, this is not moved to another location
    std::set<std::string> computeObjectList;
    RegisterComputeObjects<double, 2>(m, computeObjectList);
    RegisterComputeObjects<double, 3>(m, computeObjectList);
    RegisterComputeObjects<float, 2>(m, computeObjectList);
    RegisterComputeObjects<float, 3>(m, computeObjectList);

    RegisterComputationResult<float>(m);
    RegisterComputationResult<double>(m);

    RegisterComputations(m);

    m.attr("ComputeObjectList") = nb::cast(computeObjectList);

    // We do not bind need_other_pass or conflict found as the user has no way
    // to act on this...
    nb::enum_<Ponca::FIT_RESULT>(m, "FitResult", nb::is_arithmetic())
        .value("STABLE", Ponca::FIT_RESULT::STABLE)
        .value("UNSTABLE", Ponca::FIT_RESULT::UNSTABLE)
        .value("UNDEFINED", Ponca::FIT_RESULT::UNDEFINED)
        .export_values();
}
