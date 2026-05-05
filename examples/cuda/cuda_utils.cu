/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cuda/ponca_fit_kdtree.cu
 * \brief Utils function for using Ponca with CUDA
 * \authors Auberval Florian, Nicolas Mellado
 */

#pragma once

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>

#define CUDA_CHECK(err)                                                                                               \
    if (err != cudaSuccess)                                                                                           \
    {                                                                                                                 \
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        abort();                                                                                                      \
    }

//! [Definition KdTreeGPU]
/*! \brief A KdTree Type that can be run on the GPU
 *
 * \warning The KdTreeBase::build function cannot be used in the CUDA kernel,
 * because it still expects an STL-like container as an input.
 * This KdTree type is used to avoid the building process.
 */
template <typename DataPoint>
using KdTreeGPU = Ponca::StaticKdTreeBase<Ponca::KdTreePointerTraits<DataPoint>>;
//! [Definition KdTreeGPU]

//! [Definition KnnGraphGPU]
/*! \brief A KnnGraph Type that can be run on the GPU
 */
template <typename DataPoint>
using KnnGraphGPU = Ponca::StaticKnnGraphBase<Ponca::KnnGraphPointerTraits<DataPoint>>;
//! [Definition KnnGraphGPU]

template <typename DataPoint>
struct KdTreeRangeNeighborsFunctor
{
    static __device__ inline auto query(KdTreeGPU<DataPoint>& d_kdtree, int i, typename DataPoint::Scalar analysisScale)
        -> Ponca::KdTreeRangeIndexQuery<Ponca::KdTreePointerTraits<DataPoint>>
    {
        //! [Use KdTree.rangeNeighbors on the GPU]
        return d_kdtree.rangeNeighbors(i, analysisScale);
        //! [Use KdTree.rangeNeighbors on the GPU]
    }
};

template <typename DataPoint>
struct KdTreeKNearestNeighborsFunctor
{
    static __device__ inline auto query(KdTreeGPU<DataPoint>& d_kdtree, int i, typename DataPoint::Scalar analysisScale)
        -> Ponca::KdTreeKNearestIndexQuery<Ponca::KdTreePointerTraits<DataPoint>>
    {
        //! [Use KdTree.kNearestNeighbors on the GPU]
        return d_kdtree.kNearestNeighbors(i, d_kdtree.pointCount());
        //! [Use KdTree.kNearestNeighbors on the GPU]
    }
};

template <typename DataPoint>
struct KnnGraphRangeFunctor
{
    static __device__ inline auto query(KnnGraphGPU<DataPoint>& d_knngraph, int i,
                                        typename DataPoint::Scalar analysisScale)
        -> Ponca::KnnGraphRangeQuery<Ponca::KnnGraphPointerTraits<DataPoint>>
    {
        //! [Use KnnGraph.rangeNeighbors on the GPU]
        return d_knngraph.rangeNeighbors(i, analysisScale);
        //! [Use KnnGraph.rangeNeighbors on the GPU]
    }
};

template <typename DataPoint>
struct KnnGraphKNearestFunctor
{
    static __device__ inline auto query(KnnGraphGPU<DataPoint>& d_knngraph, int i,
                                        typename DataPoint::Scalar analysisScale)
        -> Ponca::KnnGraphKNearestQuery<Ponca::KnnGraphPointerTraits<DataPoint>>
    {
        //! [Use KnnGraph.kNearestNeighbors on the GPU]
        return d_knngraph.kNearestNeighbors(i);
        //! [Use KnnGraph.kNearestNeighbors on the GPU]
    }
};

/*! \brief Converts and uploads the internal data of the Buffers structure to the device using raw memory pointers,
 * using an intermediary binding host-device structure.
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KdTree traits structure type
 * \tparam STLBuffers A reference to a Buffer structure type that is holding STL-like containers
 * \tparam StaticBuffers A reference to a Buffer structure type that is holding memory pointers
 * \tparam IntermFunctor A function to be called right before the copy of the intermediary host structure to the device
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the
 * device
 * \param hostBuffersHoldingDevicePointers The Buffers structure on the host that references memory on the device
 * \param deviceBuffers The pointer to the Buffers structure on the device
 * \see freeBuffersOnDevice to free memory on the device with hostBuffersHoldingDevicePointers as an argument
 */
template <typename Traits, typename STLBuffers, typename StaticBuffers, typename IntermFunctor>
void deepCopyBuffersToDevice(const STLBuffers& hostBuffers, IntermFunctor midProcess, // Input
                             StaticBuffers& hostBuffersHoldingDevicePointers,
                             StaticBuffers* const deviceBuffers // Outputs
)
{
    using DataPoint = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using IndexType = typename Traits::IndexType; ///< Type used to index points into the PointContainer

    // Assign buffer sizes
    hostBuffersHoldingDevicePointers.points_size  = hostBuffers.points_size;
    hostBuffersHoldingDevicePointers.indices_size = hostBuffers.indices_size;

    // Allocate memory for the data on the device
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.points, hostBuffers.points_size * sizeof(DataPoint)));
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.indices, hostBuffers.indices_size * sizeof(IndexType)));

    // Copy the data to the device
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.points, hostBuffers.points.data(),
                          hostBuffers.points_size * sizeof(DataPoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.indices, hostBuffers.indices.data(),
                          hostBuffers.indices_size * sizeof(IndexType), cudaMemcpyHostToDevice));

    midProcess();

    // Copy host structure itself to device
    CUDA_CHECK(
        cudaMemcpy(deviceBuffers, &hostBuffersHoldingDevicePointers, sizeof(StaticBuffers), cudaMemcpyHostToDevice));
}

/*! \brief Converts and uploads the internal data of the KdTree Buffers structure to the device using raw memory
 * pointers, using an intermediary binding host-device structure.
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KdTree traits structure type
 * \tparam KdTreeBuffers A reference to a KdTree Buffer structure type that is holding STL-like containers
 * \tparam StaticKdTreeBuffers A reference to a KdTree Buffer structure type that is holding memory pointers
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the
 * device
 * \param hostBuffersHoldingDevicePointers The Buffers structure on the host that references memory on the device
 * \param deviceBuffers The pointer to the Buffers structure on the device
 * \see freeBuffersOnDevice to free memory on the device with hostBuffersHoldingDevicePointers as an argument
 */
template <typename Traits, typename KdTreeBuffers, typename StaticKdTreeBuffers>
void deepCopyKdTreeBuffersToDevice(const KdTreeBuffers& hostBuffers, // Input
                                   StaticKdTreeBuffers& hostBuffersHoldingDevicePointers,
                                   StaticKdTreeBuffers* const deviceBuffers // Outputs
)
{
    deepCopyBuffersToDevice<Traits>(
        hostBuffers,
        [&]() {
            using NodeType = typename Traits::NodeType; ///< Type of nodes used inside the KdTree
            hostBuffersHoldingDevicePointers.nodes_size = hostBuffers.nodes_size;
            CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.nodes, hostBuffers.nodes_size * sizeof(NodeType)));
            CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.nodes, hostBuffers.nodes.data(),
                                  hostBuffers.nodes_size * sizeof(NodeType), cudaMemcpyHostToDevice));
        },
        hostBuffersHoldingDevicePointers, deviceBuffers);
}

/*! \brief Converts and uploads the internal data of the KnnGraph Buffers structure to the device using raw memory
 * pointers, using an intermediary binding host-device structure.
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KnnGraph traits structure type
 * \tparam KnnGraphBuffers A reference to a KnnGraph Buffer structure type that is holding STL-like containers
 * \tparam StaticKnnGraphBuffers A reference to a KnnGraph Buffer structure type that is holding memory pointers
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the
 * device
 * \param hostBuffersHoldingDevicePointers The Buffers structure on the host that references memory on the device
 * \param deviceBuffers The pointer to the Buffers structure on the device
 * \see freeBuffersOnDevice to free memory on the device with hostBuffersHoldingDevicePointers as an argument
 */
template <typename KdTree, typename KnnGraphBuffers, typename StaticKnnGraphBuffers>
void deepCopyKnnGraphBuffersToDevice( KdTree & kdtree,
    const KnnGraphBuffers& hostBuffers, // Input
    StaticKnnGraphBuffers& hostBuffersHoldingDevicePointers,
    StaticKnnGraphBuffers* const deviceBuffers // Outputs
)
{
    using DataPoint = typename KdTree::DataPoint; ///< DataPoint given by user via Traits
    using IndexType = typename KdTree::IndexType; ///< Type used to index points into the PointContainer

    // Assign buffer sizes
    hostBuffersHoldingDevicePointers.points_size  = hostBuffers.points_size;
    hostBuffersHoldingDevicePointers.indices_size = hostBuffers.indices_size;

    // Allocate memory for the data on the device
    CUDA_CHECK(cudaMalloc((void**)&hostBuffersHoldingDevicePointers.points, hostBuffers.points_size * sizeof(DataPoint)));
    CUDA_CHECK(cudaMalloc((void**)&hostBuffersHoldingDevicePointers.indices, hostBuffers.indices_size * sizeof(IndexType)));

    // Copy the data to the device

    /////// TO FIX : Compiles only if KnnGraphPointerTraits::PointContainer is not const but "CUDA error: unspecified launch failure" at runtime
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.points, hostBuffers.points.data(), hostBuffers.points_size * sizeof(DataPoint), cudaMemcpyHostToDevice));
    //////// TO FIX : Compiles only if KnnGraphPointerTraits::PointContainer is not const but "CUDA error: unspecified launch failure" at runtime
    // CUDA_CHECK(cudaMemcpy((DataPoint*)hostBuffersHoldingDevicePointers.points, kdtree.buffers().points.data(), hostBuffers.points_size * sizeof(DataPoint), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.indices, hostBuffers.indices.data(),
                          hostBuffers.indices_size * sizeof(IndexType), cudaMemcpyHostToDevice));

    hostBuffersHoldingDevicePointers.k = hostBuffers.k;

    // Copy host structure itself to device
    CUDA_CHECK(
        cudaMemcpy(deviceBuffers, &hostBuffersHoldingDevicePointers, sizeof(StaticKnnGraphBuffers), cudaMemcpyHostToDevice));
}

/*! \brief Free the memory array internal to the Buffers on the device.
 *
 * \tparam StaticKdTreeBuffers The Buffers structure type containing the internal containers of the KdTree
 * \param hostBuffersHoldingDevicePointers The Buffers structure that references memory on the GPU
 */
template <typename StaticKdTreeBuffers>
void freeKdTreeBuffersOnDevice(const StaticKdTreeBuffers& hostBuffersHoldingDevicePointers)
{
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.points));
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.indices));
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.nodes));
}

/*! \brief Free the memory array internal to the Buffers on the device.
 *
 * \tparam StaticKnnGraphBuffers The Buffers structure type containing the internal containers of the KnnGraph
 * \param hostBuffersHoldingDevicePointers The Buffers structure that references memory on the GPU
 */
template <typename StaticKnnGraphBuffers>
void freeKnnGraphBuffersOnDevice(const StaticKnnGraphBuffers& hostBuffersHoldingDevicePointers)
{
    // CUDA_CHECK(cudaFree((Ponca::PointPositionNormal<float, 3>*)(hostBuffersHoldingDevicePointers.points))); // TO FIX : Can't be freed if of type const DataPoint*
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.points)); // TO FIX : Can't be freed if of type const DataPoint*
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.indices));
}

/*! \brief Computes the fitting process for each point of the point cloud and returns the potential and primitive
 * gradient result.
 *
 * \tparam SpatialPartitioning The Data Structure type holding the points (e.g. KdTree, KnnGraph).
 * \tparam Fit The Fit that will be computed by the Kernel.
 * \param buffers The buffers of the spatial partitioning structure.
 * \param analysisScale The radius of the neighborhood.
 * \param potentialResults As an Output, the potential results of the fit for each point of the Point Cloud.
 * \param gradientResults As an Output, the primitiveGradient results of the fit for each point of the Point Cloud.
 */
template <typename SpatialPartitioning, typename Fit, typename SpatialPartitioningQueryFunctor>
__global__ void fitPotentialAndGradientKernel(typename SpatialPartitioning::Buffers* const buffers,
                                              const typename SpatialPartitioning::DataPoint::Scalar analysisScale,
                                              typename SpatialPartitioning::DataPoint::Scalar* const potentialResults,
                                              typename SpatialPartitioning::DataPoint::Scalar* const gradientResults)
{
    using DataPoint  = typename SpatialPartitioning::DataPoint;
    using VectorType = typename DataPoint::VectorType;

    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
    // Skip when not in the point cloud
    if (i >= buffers->points_size)
        return;

    // Set up the fit
    Fit fit;
    VectorType pos = buffers->points[i].pos();
    fit.setNeighborFilter({pos, analysisScale});

    // Compute using the Spatial Partitioning data structure
    SpatialPartitioning spatialPartitioning(*buffers);
    auto rangNeighbors = SpatialPartitioningQueryFunctor::query(spatialPartitioning, i, analysisScale);
    fit.computeWithIds(rangNeighbors, buffers->points);

    // Returns NaN if not stable
    if (!fit.isStable())
    {
        potentialResults[i] = NAN;
        for (int d = 0; d < DataPoint::Dim; ++d)
        {
            gradientResults[i * DataPoint::Dim + d] = NAN;
        }
        return;
    }

    // Return the fit.potential result as an output
    potentialResults[i]   = fit.potential(pos);
    const VectorType grad = fit.primitiveGradient(pos);
    for (int d = 0; d < DataPoint::Dim; ++d)
    {
        gradientResults[i * DataPoint::Dim + d] = grad(d);
    }
}
