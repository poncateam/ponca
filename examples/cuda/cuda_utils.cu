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
 * \tparam KdTree The KdTree type
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the device
 * \return The pointer to the Buffers structure on the host that references memory on the device
 * \see freeBuffersOnDevice
 */
template <typename KdTree>
typename KdTree::Buffers deepCopyBuffersToDevice(const typename KdTree::Buffers& hostBuffers)
{
    using DataPoint      = typename KdTree::DataPoint; ///< DataPoint given by user via Traits
    using IndexType      = typename KdTree::IndexType; ///< Type used to index points into the PointContainer
    using NodeType       = typename KdTree::NodeType;  ///< Type of nodes used inside the KdTree
    using Buffers        = typename KdTree::Buffers;   ///< The Buffers structure containing the buffers stored inside the kdtree

    // The host to device copyable buffers, referencing memory on the GPU
    Buffers devBuffers{};

    // Deep copy of the internal containers
    CUDA_CHECK(cudaMalloc(&devBuffers.points,  hostBuffers.points_size  * sizeof(DataPoint)));
    CUDA_CHECK(cudaMalloc(&devBuffers.nodes,   hostBuffers.nodes_size   * sizeof(NodeType)));
    CUDA_CHECK(cudaMalloc(&devBuffers.indices, hostBuffers.indices_size * sizeof(IndexType)));

    CUDA_CHECK(cudaMemcpy(devBuffers.points,  hostBuffers.points,
               hostBuffers.points_size * sizeof(DataPoint),
               cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(devBuffers.nodes,   hostBuffers.nodes,
               hostBuffers.nodes_size * sizeof(NodeType),
               cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(devBuffers.indices, hostBuffers.indices,
               hostBuffers.indices_size * sizeof(IndexType),
               cudaMemcpyHostToDevice));

    devBuffers.points_size  = hostBuffers.points_size;
    devBuffers.nodes_size   = hostBuffers.nodes_size;
    devBuffers.indices_size = hostBuffers.indices_size;

    return devBuffers;
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
 * \warning The build function cannot be used in the CUDA kernel,
 * because it still expects an STL-like container as an input.
 * This KdTree must be instanciated with the `KdTreeBase::KdTreeBase` constructor on the GPU.
 */
template <typename DataPoint>
using KdTreeGPU = Ponca::KdTreeDenseBase<Ponca::KdTreePointerTraits<DataPoint>>;
