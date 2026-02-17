#pragma once

#define CUDA_CHECK(err) \
    if (err != cudaSuccess) { \
    std::cerr << "CUDA error: " << cudaGetErrorString(err) \
    << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    abort(); \
}

/*! \brief Kernel that binds each point buffer to the interlaced buffer of positions and normals in the Device.
 *
 * \tparam DataPoint The DataPoint type.
 * \tparam Scalar A scalar type (float, double...)
 * \param points The point cloud (needs to be bind to the interlacedArray on the device).
 * \param interlacedArray The interlaced buffer of positions and normals.
 * \param nbPoints The total number of points in the point cloud.
 */
template<typename DataPoint, typename Scalar>
__global__ void bindPointsKernel(
    DataPoint* points,
    Scalar* interlacedArray,
    const int nbPoints
) {
    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

    // Skip when not in the point cloud
    if (i >= nbPoints) return;

    points[i].bind(interlacedArray);
}

/*!
 * \brief The default traits type used by the kd-tree.
 *
 * \see KdTreeCustomizableNode Helper class to modify Inner/Leaf nodes without redefining a Trait class
 *
 * \tparam _NodeType Type used to store nodes, set by default to #KdTreeDefaultNode
 */
template <typename _DataPoint,
        template <typename /*Index*/,
                  typename /*NodeIndex*/,
                  typename /*DataPoint*/,
                  typename /*LeafSize*/> typename _NodeType = Ponca::KdTreeDefaultNode>
struct KdTreeGPUTraits
{
    enum
    {
        /*!
         * \brief A compile-time constant specifying the maximum depth of the kd-tree.
         */
        MAX_DEPTH = 32,
    };

    /*!
     * \brief The type used to store point data.
     *
     * Must provide `Scalar` and `VectorType` typedefs.
     *
     * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
     * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
     * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
     */
    using DataPoint    = _DataPoint;
    using IndexType    = int;
    using LeafSizeType = unsigned short;

    // Containers
    using PointContainer = DataPoint*;
    using IndexContainer = IndexType*;

    // Nodes
    using NodeIndexType = std::size_t;
    using NodeType      = _NodeType<IndexType, NodeIndexType, DataPoint, LeafSizeType>;
    using NodeContainer = NodeType*;

    /*!
     * \brief Converts the STL-Like input container to a memory array to be used internally.
     *
     * \see KdTreeBase
     */
    template <typename InternalContainer, typename InputContainer>
    [[nodiscard]] static PONCA_MULTIARCH_HOST inline InternalContainer& toInternalContainer ( InputContainer & input)
    {
        return input.data();
    }
};

template <typename DataPoint>
using KdTreeGPU = Ponca::KdTreeBase<KdTreeGPUTraits<DataPoint>>;
