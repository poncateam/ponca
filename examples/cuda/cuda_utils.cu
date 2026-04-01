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

#define CUDA_CHECK(err) \
    if (err != cudaSuccess) { \
    std::cerr << "CUDA error: " << cudaGetErrorString(err) \
    << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    abort(); \
}

/*! \brief Uploads the internal data of the Buffers structure to the device,
 * using an intermediary binding host-device structure.
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KdTree traits structure type
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the device
 * \param hostBuffersHoldingDevicePointers The Buffers structure on the host that references memory on the device
 * \param deviceBuffers The pointer to the Buffers structure on the device
 * \see freeBuffersOnDevice to free memory on the device with hostBuffersHoldingDevicePointers as an argument
 */
template <typename Traits, bool SendNode = true, typename KdTreeDenseBuffers, typename StaticKdTreeBuffers>
void deepCopyBuffersToDevice(const KdTreeDenseBuffers& hostBuffers, StaticKdTreeBuffers& hostBuffersHoldingDevicePointers, StaticKdTreeBuffers* const deviceBuffers)
{
    using DataPoint      = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using IndexType      = typename Traits::IndexType; ///< Type used to index points into the PointContainer

    // Assign buffer sizes
    hostBuffersHoldingDevicePointers.points_size  = hostBuffers.points_size;
    hostBuffersHoldingDevicePointers.indices_size = hostBuffers.indices_size;

    // Allocate memory for the data on the device
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.points,  hostBuffers.points_size  * sizeof(DataPoint)));
    CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.indices, hostBuffers.indices_size * sizeof(IndexType)));

    // Copy the data to the device
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.points, hostBuffers.points.data(),
        hostBuffers.points_size * sizeof(DataPoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.indices, hostBuffers.indices.data(),
        hostBuffers.indices_size * sizeof(IndexType), cudaMemcpyHostToDevice));
    if constexpr (SendNode) {
        using NodeType       = typename Traits::NodeType;  ///< Type of nodes used inside the KdTree
        hostBuffersHoldingDevicePointers.nodes_size   = hostBuffers.nodes_size;
        CUDA_CHECK(cudaMalloc(&hostBuffersHoldingDevicePointers.nodes,   hostBuffers.nodes_size   * sizeof(NodeType)));
        CUDA_CHECK(cudaMemcpy(hostBuffersHoldingDevicePointers.nodes, hostBuffers.nodes.data(),
            hostBuffers.nodes_size * sizeof(NodeType), cudaMemcpyHostToDevice));
    }

    // Copy host structure itself to device
    CUDA_CHECK(cudaMemcpy(deviceBuffers, &hostBuffersHoldingDevicePointers,
        sizeof(StaticKdTreeBuffers), cudaMemcpyHostToDevice)
    );
}

/*! \brief Free the memory array internal to the KdTree Buffers on the device.
 *
 * \tparam Buffers The Buffers structure type containing the internal containers of the KdTree
 * \param hostBuffersHoldingDevicePointers The Buffers structure that references memory on the GPU
 */
template <bool SendNode = true, typename Buffers>
void freeBuffersOnDevice(const Buffers& hostBuffersHoldingDevicePointers)
{
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.points));
    CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.indices));
    if constexpr (SendNode)
        CUDA_CHECK(cudaFree(hostBuffersHoldingDevicePointers.nodes));
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
