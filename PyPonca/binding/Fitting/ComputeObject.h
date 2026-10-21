/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "../Common/pypoint.h"
#include "../Common/DeviceUtils.h"

/**
 * \brief Describes a computation to be performed
 */
template <typename Scalar>
struct ComputationDescriptor
{
    /**
     * Id of the computation (likely bound to an enum, but we this is more general)
     */
    size_t id;

    /**
     * Expected output dimension. The first argument will be overriden by the amount of
     * filters required.
     */
    std::vector<size_t> outputDims;

    /**
     * Input data
     */
    nb::ndarray<Scalar> inputData;
};

/**
 * \brief Describes the result of a computation
 *
 * For now, this structure only contain the result array.
 * Additionnal information can be infered by the user through
 * the order in which this will be returned.
 */

template <typename Scalar>
struct ComputationResults
{
    nb::ndarray<int, nb::array_api> fitStatus;
    std::vector<nb::ndarray<Scalar, nb::array_api>> results;
};

/**
 * \brief Binding to a Ponca::ComputeObject
 */
template <typename _CO>
struct PyComputeObject
{
public:
    using ComputeObject = _CO;
    using Point         = typename ComputeObject::DataPoint;
    using Scalar        = typename Point::Scalar;
    using VectorType    = typename Point::VectorType;

    /**
     * \brief Default constructor
     */
    PyComputeObject() { m_hasFilter = false; }

    /**
     * \brief Sets the neighborFilter informations
     *
     * Nothing will be set for now, informations are stored for further usage
     *
     * \param _pos The position array
     * \param _t The scalar (radius) array
     */
    void setNeighborFilter(const PyVectorArray<Point>& _pos, const PyScalarArray<Point>& _t,
                           const PyVectorArray<Point>& _normals)
    {
        if (_pos.device_type() != _t.device_type())
            throw std::invalid_argument("Incompatible devices between pos and radii.");

        if (_pos.ndim() != 2)
            throw std::invalid_argument("Only 2D array are supported for filter centers");

        if (_t.ndim() != 1)
            throw std::invalid_argument("Only scalar array are supported for filter radii");

        if (_t.shape(0) != _pos.shape(0) && _t.shape(0) != 1)
            throw std::invalid_argument("Incompatible shapes between filter pos and radii");

        if constexpr (RequiresFilterNormals())
        {
            if (_normals.data() == nullptr)
                throw std::invalid_argument("This compute object requires normals.");

            if (_normals.device_type() != _pos.device_type())
                throw std::invalid_argument("Incompatible devices between pos and normals.");

            if (_normals.ndim() != 2)
                throw std::invalid_argument("Only 2D array are supported for filter normals");

            if (_normals.shape(0) != _pos.shape(0))
                throw std::invalid_argument("Shape mismatch between pos and normals");

            if (_normals.shape(1) != _pos.shape(1))
                throw std::invalid_argument("Shape mismatch between pos and normals");
        }

        m_filterCenters = _pos;
        m_filterRadius  = PyScalarArray<Point>(broadcastArrayTo<Scalar, 1>(_t, {_pos.shape(0)}));
        m_filterNormals = _normals;

        m_hasFilter = true;
    }

    /**
     * \brief Adds a computation to be performed
     *
     * \param descriptor The description
     */
    void addComputation(ComputationDescriptor<Scalar>&& descriptor) { m_descriptors.push_back(std::move(descriptor)); }

    /**
     * \brief Perform the computation
     *
     * The compute functor first defines how the point cloud is used and fed to the computation.
     * This allows this function to remain agnostic to the cloud representation.
     * The expected signature is:
     * - Co& co: the compute cobject. Do not forget the reference !!!
     * - const Cloud& cloud: The data representation of the point cloud
     * - unsigned int i: The index of the filter
     * - Vector loc: The filter location
     * - Scalar t: The filter radius
     *
     * A function is given in order to extract results once the main compute is performed.
     * The expected signature is:
     *  - unsigned int id: the id of the computation
     *  - CO& co: the fitted compute object
     *  - unsigned int i: the current center index
     *  - const ArrayView& input: input array for the computation
     *  - ArrayView& output: output array for the computation
     *
     * \param cloud The Point cloud
     * \param f Function to extract result
     */
    template <typename PointCloud, typename ComputeFunctor, typename ExtractFunc>
    ComputationResults<Scalar> compute(const PointCloud& cloud, ComputeFunctor&& cf, ExtractFunc&& ef)
    {
        // Verify integrity of data
        if (m_descriptors.size() == 0)
            throw std::runtime_error("No computation to be performed");

        if (!m_hasFilter)
            throw std::runtime_error("No filter set for this computation");

        if (cloud.device_type() != m_filterCenters.device_type())
            throw std::runtime_error("Device mismatch between point cloud and filters");

        for (size_t i = 0; i < m_descriptors.size(); ++i)
            if (cloud.device_type() != m_descriptors[i].inputData.device_type() &&
                m_descriptors[i].inputData.data() != nullptr)
                throw std::runtime_error("Device mismatch between computation input " + std::to_string(i) +
                                         " and point cloud.");

        // Fill in the first dimension of outputs
        for (size_t i = 0; i < m_descriptors.size(); ++i)
            if (m_descriptors[i].outputDims.size() > 0)
                m_descriptors[i].outputDims[0] = m_filterCenters.shape(0);

        switch (cloud.device_type())
        {
        case nb::device::cpu::value:
            return computeCPU(cloud, std::forward<ComputeFunctor>(cf), std::forward<ExtractFunc>(ef));
#ifdef __CUDACC__
//         case nb::device::cuda::value:
//         return computeCUDA(cloud, std::forward<ComputeFunctor>(cf), std::forward<ExtractFunc>(ef));
#endif
        default:
            throw std::runtime_error("Unsupported device: " + std::to_string(cloud.device_type()));
        }

        return {};
    }

    static constexpr bool RequiresFilterNormals()
    {
        using NF = typename ComputeObject::NeighborFilter;
        return requires { NF(VectorType{}, Scalar{}, VectorType{}); };
    }

private:
    /**
     * \brief Creates a filter specialized for
     */
    auto MakeFilter(size_t i) const
    {
        using NF                 = typename ComputeObject::NeighborFilter;
        const auto queryLocation = PyVectorArrayIndex<Point>(m_filterCenters, i);
        const auto queryRadius   = PyScalarArrayIndex<Point>(m_filterRadius, i);

        if constexpr (RequiresFilterNormals())
        {
            const auto queryNormal = PyVectorArrayIndex<Point>(m_filterNormals, i);
            return NF(queryLocation, queryRadius, queryNormal);
        }
        else
            return NF(queryLocation, queryRadius);
    }

    /**
     * \brief Perform the computation on the CPU
     *
     * The compute functor first defines how the point cloud is used and fed to the computation.
     * This allows this function to remain agnostic to the cloud representation.
     * The expected signature is:
     * - Co& co: the compute cobject. Do not forget the reference !!!
     * - const Cloud& cloud: The data representation of the point cloud
     * - unsigned int i: The index of the filter
     * - Co::Filter nf: The neighbor filter
     *
     * A function is given in order to extract results once the main compute is performed.
     * The expected signature is:
     *  - unsigned int id: the id of the computation
     *  - CO& co: the fitted compute object
     *  - unsigned int i: the current center index
     *  - const ArrayView& input: input array for the computation
     *  - ArrayView& output: output array for the computation
     *
     * \param cloud The Point cloud representation
     * \param cf The compute functor. This function must be specialized for each representation of the point cloud (ie.
     * array vs kdtree)
     * \param ef Function to extract result
     */
    template <typename PointCloud, typename ComputeFunctor, typename ExtractFunc>
    ComputationResults<Scalar> computeCPU(const PointCloud& cloud, ComputeFunctor&& cf, ExtractFunc&& ef)
    {
        ComputationResults<Scalar> results;

        DeviceArrayView<int> statusView;
        std::vector<DeviceArrayView<Scalar>> inputsView;
        std::vector<DeviceArrayView<Scalar>> resultView;

        results.fitStatus = CreateDeviceArray<int>({m_filterCenters.shape(0)}, cloud.device_type());
        statusView        = nb::ndarray<int>(results.fitStatus);
        for (size_t i = 0; i < m_descriptors.size(); ++i)
        {
            results.results.emplace_back(CreateDeviceArray<Scalar>(m_descriptors[i].outputDims, cloud.device_type()));

            inputsView.emplace_back(m_descriptors[i].inputData);
            resultView.emplace_back(nb::ndarray<Scalar>(results.results.back()));
        }

        for (unsigned int i = 0; i < m_filterCenters.shape(0); ++i)
        {
            ComputeObject computeObject;

            computeObject.init();
            computeObject.setNeighborFilter(MakeFilter(i));
            statusView.data[i] = static_cast<int>(cf(computeObject, cloud, i));

            // We do not check for result. The cost of running this loop is low and avoid
            // unnecessary if.
            for (size_t j = 0; j < m_descriptors.size(); ++j)
                ef(m_descriptors[j].id, computeObject, i, inputsView[j], resultView[j]);
        }

        return results;
    }

private:
    /**
     * Stores if setNeighborFilter was called
     */
    bool m_hasFilter;

    /**
     * Computation to be performed
     */
    std::vector<ComputationDescriptor<Scalar>> m_descriptors;

    /**
     * Filter centers
     */
    PyVectorArray<Point> m_filterCenters;
    /**
     * Filers radiuses
     */
    PyScalarArray<Point> m_filterRadius;
    /**
     * Filters normals
     */
    PyVectorArray<Point> m_filterNormals;
};

