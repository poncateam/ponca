#pragma once

#define CUDA_CHECK(err) \
    if (err != cudaSuccess) { \
    std::cerr << "CUDA error: " << cudaGetErrorString(err) \
    << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    abort(); \
}

/*! \brief A KdTree Type that can be run on the GPU
 *
 * \warning The build function cannot be used in the CUDA kernel,
 * because it still expects an STL-like container as an input.
 * This KdTree must be instanciated with the `KdTreeBase::KdTreeBase` constructor on the GPU.
 */
template <typename DataPoint>
using KdTreeGPU = Ponca::KdTreeBase<Ponca::KdTreePointerTraits<DataPoint>>;
