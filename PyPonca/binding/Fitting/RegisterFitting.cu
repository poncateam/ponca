/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <algorithm>
#include <string>

#include "FittingList.h"
#include "../SpatialPartitioning/pykdtree.h"

// This file is the main binding code for the Fitting module.
//
// The philosophy is to avoid lambdas as much as possible and rely
// on plain function instead. Lambdas may or may not be supported by
// some device code. In CUDA for instance, only one layer of lambda
// is supported.

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
 * 
 * \tparam Co The compute obejct to perform the compute method on
 * \tparam Cloud Pointcloud type
 * 
 * \param co The compute object
 * \param i Index into the filter list
 * \param loc Filter center
 * \param rad Filter radius
 */
template <typename Co, typename Cloud>
void PerformRawCloudComputation(Co& co, const Cloud& cloud, unsigned int i, typename Cloud::Point::VectorType loc,
                                typename Cloud::Point::Scalar rad)
{
    co.compute(cloud.begin(), cloud.end());
}

/**
 * \brief Perform the compute method on a KDTree
 * 
 * \tparam Co The compute obejct to perform the compute method on
 * \tparam Cloud Pointcloud type
 * 
 * \param co The compute object
 * \param i Index into the filter list
 * \param loc Filter center
 * \param rad Filter radius
 */
template <typename Co, typename PyKDTree>
void PerformKDTreeComputation(Co& co, const PyKDTree& kdtree, unsigned int i, typename PyKDTree::Point::VectorType loc,
                              typename PyKDTree::Point::Scalar rad)
{
    kdtree.Run([&](const auto& tree) {
        auto neighbors = tree.rangeNeighbors(loc, rad);

        std::vector<int> indices;
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(indices));
        co.computeWithIds(indices, tree.points());
    });
}

/**
 * \brief Performs the computation
 *
 * \param self The instance of the object
 * \param cloud The point cloud
 */
template <typename PyCo, typename Cloud>
auto computeRawPointCloud(PyCo& self, const Cloud& cloud)
{
    using CO = typename PyCo::ComputeObject;
    return self.compute(cloud, &PerformRawCloudComputation<CO, Cloud>, &ExtractComputation<CO>);
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
    using CO = typename PyCo::ComputeObject;
    return self.compute(cloud, &PerformKDTreeComputation<CO, PyKDTree>, &ExtractComputation<CO>);
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

    // Compute mangling informations
    const std::string mangledName = PointCloud::PointName + _NF::name;

    // General properties
    Factory::foreach ([&](const auto& x) {
        using T                   = decltype(x.object);
        using PyCo                = PyComputeObject<T>;
        
        std::string coname = x.name;
        coname.erase(std::remove(coname.begin(), coname.end(), ' '), coname.end());

        const std::string newname = coname + mangledName;
        auto pyco = nb::class_<PyCo>(m, newname.c_str());
        pyco.def(nb::init<>());
        pyco.def("setNeighborFilter", &PyCo::setNeighborFilter);
        pyco.def("addComputation", &addComputation<PyCo>, nb::arg("id"), nb::arg("data").none());
        pyco.def("compute", &computeRawPointCloud<PyCo, PointCloud>);
        pyco.def("compute", &computeKDTree<PyCo, PyKDTree<PointCloud>>);
        
        list.insert(coname);
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
    const std::string computationResultName = "ComputationResultlt" + MangleType<Scalar>();

    auto result = nb::class_<ComputationResult<Scalar>>(m, computationResultName.c_str());
    result.def_rw("data", &ComputationResult<Scalar>::resultData, nb::rv_policy::reference);
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
}
