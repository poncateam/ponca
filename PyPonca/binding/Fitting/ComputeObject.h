/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "../Common/pypoint.h"
#include "../Common/DeviceUtils.h"

#include <optional>

#ifdef PYPONCA_WITH_OMP
#    include <omp.h>
#endif

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
    /**
     * Status for each fit
     */
    nb::ndarray<int, nb::array_api> fitStatus;

    /**
     * Result for each computation, in the order first provided
     */
    std::vector<nb::ndarray<Scalar, nb::array_api>> results;
};

/**
 * \brief Device-agnostic views to data used during a computation
 */
template <typename Scalar>
struct ComputationViews
{
    // Number of views describing the filters
    static constexpr size_t FilterViewCount = 3;

    /**
     * \brief Ctor
     *
     * \param views Views over filter centers, radii and normals, followed by the input of each
     * computation, then by the output of each computation (FilterViewCount + 2 * computationCount views)
     * \param ids Id of each computation
     * \param computationCount Number of computations
     * \param status Fit status of each filter
     */
    PONCA_MULTIARCH ComputationViews(DeviceArrayView<Scalar>* views, const size_t* ids, size_t computationCount,
                                     int* status)
        : filters(views), computationCount(computationCount), ids(ids), inputs(views + FilterViewCount),
          outputs(views + FilterViewCount + computationCount), status(status)
    {
    }

    /**
     * Filters data (in this order: center, radius, normals)
     */
    const DeviceArrayView<Scalar>* filters;

    /**
     * Number of computations = size of the arrays
     */
    size_t computationCount;

    /**
     * Computation ids
     */
    const size_t* ids;

    /**
     * Input data
     */
    const DeviceArrayView<Scalar>* inputs;

    /**
     * Output data
     */
    DeviceArrayView<Scalar>* outputs;

    /**
     * Output status
     */
    int* status;
};

#ifdef __CUDACC__
/**
 * \brief Entry point for cuda compute.
 *
 * \see PyComputeObject::ComputeFilters
 */
template <typename PyCo, typename Cloud, typename ComputeFunctor, typename ExtractFunc, typename Scalar>
__global__ void ComputeKernel(const Cloud cloud, const ComputeFunctor cf, const ExtractFunc ef,
                              const ComputationViews<Scalar> views, size_t filterCount)
{
    const size_t first = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t step  = static_cast<size_t>(gridDim.x) * blockDim.x;
    PyCo::Compute(first, step, cloud, cf, ef, views, filterCount);
}
#endif

/**
 * \brief Binding to a Ponca::ComputeObject
 *
 * \tparam _CO The compute object
 * \tparam _DeviceSupport Whether the compute object can be used within device code. When false,
 * computations on devices raise an error, and no kernel is instantiated.
 */
template <typename _CO, bool _DeviceSupport>
struct PyComputeObject
{
public:
    using ComputeObject = _CO;
    using Point         = typename ComputeObject::DataPoint;
    using Scalar        = typename Point::Scalar;
    using VectorType    = typename Point::VectorType;

    static constexpr bool DeviceSupport = _DeviceSupport;

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
     *
     * The expected signature is:
     * - Co& co: the compute cobject. Do not forget the reference !!!
     * - const Cloud& cloud: The view of the point cloud
     * - unsigned int i: The index of the filter
     *
     * A function is given in order to extract results once the main compute is performed.
     * The expected signature is:
     *  - unsigned int id: the id of the computation
     *  - CO& co: the fitted compute object
     *  - unsigned int i: the current center index
     *  - const ArrayView& input: input array for the computation
     *  - ArrayView& output: output array for the computation
     *
     * Both functors must be callable from device code (see PONCA_MULTIARCH).
     *
     * \param cloud The Point cloud
     * \param cf The compute functor
     * \param ef The extract functor
     */
    template <typename PointCloud, typename ComputeFunctor, typename ExtractFunc>
    ComputationResults<Scalar> compute(const PointCloud& cloud, ComputeFunctor&& cf, ExtractFunc&& ef)
    {
        const int device   = cloud.device_type();
        const int deviceId = cloud.device_id();

        // Verify integrity of data
        if (m_descriptors.size() == 0)
            throw std::runtime_error("No computation to be performed");

        if (!m_hasFilter)
            throw std::runtime_error("No filter set for this computation");

        if (device != m_filterCenters.device_type() || deviceId != m_filterCenters.device_id())
            throw std::runtime_error("Device mismatch between point cloud and filters");

        for (size_t i = 0; i < m_descriptors.size(); ++i)
            if ((device != m_descriptors[i].inputData.device_type() ||
                 deviceId != m_descriptors[i].inputData.device_id()) &&
                m_descriptors[i].inputData.data() != nullptr)
                throw std::runtime_error("Device mismatch between computation input " + std::to_string(i) +
                                         " and point cloud.");

        switch (device)
        {
        case nb::device::cpu::value:
            break;
#ifdef __CUDACC__
        case nb::device::cuda::value:
            if (!DeviceSupport)
                throw std::runtime_error(
                    "This compute object is not available on CUDA devices, use CPU arrays instead.");
            break;
#endif
        default:
            throw std::runtime_error("Unsupported device: " + std::to_string(device));
        }

        // Fill in the first dimension of outputs
        const size_t filterCount = m_filterCenters.shape(0);
        for (size_t i = 0; i < m_descriptors.size(); ++i)
            if (m_descriptors[i].outputDims.size() > 0)
                m_descriptors[i].outputDims[0] = filterCount;

#ifdef __CUDACC__
        // Allocations are performed on the current device
        std::optional<CudaDeviceGuard> guard;
        if (device == nb::device::cuda::value)
            guard.emplace(deviceId);
#endif

        ComputationResults<Scalar> results;
        results.fitStatus = CreateDeviceArray<int>({filterCount}, device, deviceId);

        std::vector<size_t> ids;
        for (size_t i = 0; i < m_descriptors.size(); ++i)
        {
            ids.push_back(m_descriptors[i].id);
            results.results.emplace_back(CreateDeviceArray<Scalar>(m_descriptors[i].outputDims, device, deviceId));
        }

        // Create views for compatibility with devices.
        std::vector<DeviceArrayView<Scalar>> views;
        views.reserve(ComputationViews<Scalar>::FilterViewCount + 2 * m_descriptors.size());
        views.emplace_back(nb::ndarray<Scalar>(m_filterCenters));
        views.emplace_back(nb::ndarray<Scalar>(m_filterRadius));
        views.emplace_back(nb::ndarray<Scalar>(m_filterNormals));

        for (size_t i = 0; i < m_descriptors.size(); ++i)
            views.emplace_back(m_descriptors[i].inputData);

        for (size_t i = 0; i < m_descriptors.size(); ++i)
            views.emplace_back(nb::ndarray<Scalar>(results.results[i]));

        // Computations access bytewise copies of the views and ids, stored on their device
        const auto viewsArray = CopyToDevice(views, device, deviceId);
        const auto idsArray   = CopyToDevice(ids, device, deviceId);
        const ComputationViews<Scalar> computation(reinterpret_cast<DeviceArrayView<Scalar>*>(viewsArray.data()),
                                                   reinterpret_cast<const size_t*>(idsArray.data()), ids.size(),
                                                   results.fitStatus.data());

        // Only the distribution of the filters depends on the device
        const auto cloudView = cloud.view();
        switch (device)
        {
        case nb::device::cpu::value:
            // As on the device, each thread processes a part of the filters
#pragma omp parallel
        {
#ifdef PYPONCA_WITH_OMP
            const size_t first = omp_get_thread_num();
            const size_t step  = omp_get_num_threads();
#else
            const size_t first = 0;
            const size_t step  = 1;
#endif
            Compute(first, step, cloudView, cf, ef, computation, filterCount);
        }
        break;
#ifdef __CUDACC__
        case nb::device::cuda::value:
            if constexpr (DeviceSupport)
            {
                if (filterCount == 0)
                    break;

                constexpr unsigned int blockSize = 128;
                const unsigned int gridSize      = static_cast<unsigned int>((filterCount + blockSize - 1) / blockSize);

                ComputeKernel<PyComputeObject><<<gridSize, blockSize>>>(cloudView, cf, ef, computation, filterCount);

                CheckCuda(cudaGetLastError(), "Unable to launch the computation kernel");
                CheckCuda(cudaDeviceSynchronize(), "Computation kernel failed");
            }
            break;
#endif
        }

        return results;
    }

    static constexpr bool RequiresFilterNormals()
    {
        using NF = typename ComputeObject::NeighborFilter;
        return requires { NF(VectorType{}, Scalar{}, VectorType{}); };
    }

    /**
     * \brief Fits the filters i = first + k * step (i < filterCount) and extracts the results
     * of every computation
     *
     * This is the computation loop, whatever the device: each thread (OpenMP thread on the CPU,
     * thread of ComputeKernel on the device) processes a part of the filters.
     *
     * \param first Index of the first filter to process
     * \param step Distance between two filters to process
     * \param cloud View of the point cloud (see compute)
     * \param cf The compute functor (see compute)
     * \param ef The extract functor (see compute)
     * \param views Arrays of the computation
     * \param filterCount Number of filters
     */
    template <typename Cloud, typename ComputeFunctor, typename ExtractFunc>
    static PONCA_MULTIARCH void Compute(size_t first, size_t step, const Cloud& cloud, const ComputeFunctor& cf,
                                        const ExtractFunc& ef, const ComputationViews<Scalar>& views,
                                        size_t filterCount)
    {
        using NF = typename ComputeObject::NeighborFilter;

        for (size_t i = first; i < filterCount; i += step)
        {
            const VectorType queryLocation = views.filters[0].template GetVector<Point>(i);
            const Scalar queryRadius       = *(views.filters[1].data + i * views.filters[1].stride[0]);

            ComputeObject computeObject;
            computeObject.init();
            if constexpr (RequiresFilterNormals())
                computeObject.setNeighborFilter(
                    NF(queryLocation, queryRadius, views.filters[2].template GetVector<Point>(i)));
            else
                computeObject.setNeighborFilter(NF(queryLocation, queryRadius));

            views.status[i] = static_cast<int>(cf(computeObject, cloud, i));

            // We do not check for result. The cost of running this loop is low and avoid
            // unnecessary if.
            for (size_t j = 0; j < views.computationCount; ++j)
                ef(views.ids[j], computeObject, i, views.inputs[j], views.outputs[j]);
        }
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
