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

/*! \brief A KdTree Type that can be run on the GPU
 *
 * \warning The build function cannot be used in the CUDA kernel,
 * because it still expects an STL-like container as an input.
 * This KdTree must be instanciated with the `KdTreeBase::KdTreeBase` constructor on the GPU.
 */
template <typename DataPoint>
using KdTreeGPU = Ponca::KdTreeBase<Ponca::KdTreePointerTraits<DataPoint>>;
