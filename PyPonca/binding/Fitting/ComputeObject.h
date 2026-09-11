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
struct ComputationResult
{
    nb::ndarray<Scalar, nb::array_api> resultData;
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
    void setNeighborFilter(const PyVectorArray<Point>& _pos, const PyScalarArray<Point>& _t)
    {
        if (_pos.device_type() != _t.device_type())
            throw std::runtime_error("Incompatible devices between pos and t.");

        if (_pos.ndim() != 2)
            throw std::runtime_error("Only 2D array are supported for filter centers");

        if (_t.ndim() != 1)
            throw std::runtime_error("Only scalar array are supported for filter radii");

        m_filterCenters = _pos;
        m_filterRadius  = PyScalarArray<Point>(broadcastArrayTo<Scalar, 1>(_t, {_pos.shape(0)}));

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
    std::vector<ComputationResult<Scalar>> compute(const PointCloud& cloud, ComputeFunctor&& cf, ExtractFunc&& ef)
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

private:
    /**
     * \brief Perform the computation on the CPU
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
     * \param cloud The Point cloud representation
     * \param cf The compute functor. This function must be specialized for each representation of the point cloud (ie.
     * array vs kdtree)
     * \param ef Function to extract result
     */
    template <typename PointCloud, typename ComputeFunctor, typename ExtractFunc>
    std::vector<ComputationResult<Scalar>> computeCPU(const PointCloud& cloud, ComputeFunctor&& cf, ExtractFunc&& ef)
    {
        std::vector<ComputationResult<Scalar>> result;
        std::vector<DeviceArrayView<Scalar>> inputsView;
        std::vector<DeviceArrayView<Scalar>> resultView;

        // Prepare array views for computation
        for (size_t i = 0; i < m_descriptors.size(); ++i)
        {
            ComputationResult<Scalar> computationResult;
            computationResult.resultData = CreateDeviceArray<Scalar>(m_descriptors[i].outputDims, cloud.device_type());

            inputsView.emplace_back(m_descriptors[i].inputData);
            resultView.emplace_back(nb::ndarray<Scalar>(computationResult.resultData));
                
            result.push_back(computationResult);
        }
        
        for (unsigned int i = 0; i < m_filterCenters.shape(0); ++i)
        {
            ComputeObject computeObject;

            const auto queryLocation = PyVectorArrayIndex<Point>(m_filterCenters, i);
            const auto queryRadius   = PyScalarArrayIndex<Point>(m_filterRadius, i);

            computeObject.init();
            computeObject.setNeighborFilter({queryLocation, queryRadius});
            cf(computeObject, cloud, i, queryLocation, queryRadius);

            for (size_t j = 0; j < m_descriptors.size(); ++j)
                ef(m_descriptors[j].id, computeObject, i, inputsView[j], resultView[j]);
        }

        return result;
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
};

