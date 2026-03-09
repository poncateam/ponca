#pragma once

#define CUDA_CHECK(err) \
    if (err != cudaSuccess) { \
    std::cerr << "CUDA error: " << cudaGetErrorString(err) \
    << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    abort(); \
}

/*! \brief Uploads the internal containers of the Buffers structure to the device
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KdTree traits structure type
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the device
 * \param hostBuffersHoldingDevicePointers The pointer to the Buffers structure on the host that references memory on the device
 * \param deviceBuffers The Buffers structure on the device
 * \see freeBuffersOnDevice
 */
template <typename Traits, typename KdTreeDenseBuffers, typename StaticKdTreeBuffers>
void deepCopyBuffersToDevice(const KdTreeDenseBuffers& hostBuffers, StaticKdTreeBuffers& hostBuffersHoldingDevicePointers, StaticKdTreeBuffers* deviceBuffers)
{
    using DataPoint      = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using IndexType      = typename Traits::IndexType; ///< Type used to index points into the PointContainer
    using NodeType       = typename Traits::NodeType;  ///< Type of nodes used inside the KdTree

    // Assign buffer sizes
    hostBuffersHoldingDevicePointers.points_size  = hostBuffers.points_size;
    hostBuffersHoldingDevicePointers.nodes_size   = hostBuffers.nodes_size;
    hostBuffersHoldingDevicePointers.indices_size = hostBuffers.indices_size;

    // Allocate memory for the data on the device
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.points,  hostBuffers.points_size  * sizeof(DataPoint)));
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.nodes,   hostBuffers.nodes_size   * sizeof(NodeType)));
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.indices, hostBuffers.indices_size * sizeof(IndexType)));

    // Copy the data to the device
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.points, hostBuffers.points.data(),
        hostBuffers.points_size * sizeof(DataPoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.nodes, hostBuffers.nodes.data(),
        hostBuffers.nodes_size * sizeof(NodeType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.indices, hostBuffers.indices.data(),
        hostBuffers.indices_size * sizeof(IndexType), cudaMemcpyHostToDevice));

    // Copy host structure itself to device
    CUDA_CHECK(cudaMemcpy(deviceBuffers, &hostBuffersHoldingDevicePointers,
        sizeof(StaticKdTreeBuffers), cudaMemcpyHostToDevice)
    );
}

/*! \brief Free the memory array internal to the KdTree Buffers on the device
 *
 * \tparam Buffers The Buffers structure type containing the internal containers of the KdTree
 * \param hostBuffers The Buffers structure that references memory on the GPU
 */
template <typename Buffers>
void freeBuffersOnDevice(const Buffers& hostBuffers)
{
    CUDA_CHECK(cudaFree(hostBuffers.points));
    CUDA_CHECK(cudaFree(hostBuffers.nodes));
    CUDA_CHECK(cudaFree(hostBuffers.indices));
}

/*! \brief A KdTree Type that can be run on the GPU
 *
 * \warning The KdTreeBase::build function cannot be used in the CUDA kernel,
 * because it still expects an STL-like container as an input.
 * This KdTree type is used to avoid the building process.
 */
template <typename DataPoint>
using KdTreeGPU = Ponca::StaticKdTreeBase<Ponca::KdTreePointerTraits<DataPoint>>;
