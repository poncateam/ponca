/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cuda/ponca_fit_kdtree.cu
 * \brief Example that uses the Fitting and SpatialPartitioning module with Sycl
 * \authors Bastien DOIGNIES
 */

#pragma once

#include <Ponca/Ponca>
#include <sycl/sycl.hpp>

/*! \brief Uploads the internal data of the Buffers structure to the device,
 * using an intermediary binding host-device structure.
 *
 * \warning The buffers needs to be freed after use, with the \ref freeBuffersOnDevice function
 *
 * \tparam Traits The KdTree traits structure type
 * \param queue Sycl queue
 * \param hostBuffers The Buffers structure holding the internal storage of the KdTree that we want to upload to the
 * device
 * \param hostBuffersHoldingDevicePointers The Buffers structure on the host that references memory on the device
 * \see freeBuffersOnDevice to free memory on the device with hostBuffersHoldingDevicePointers as an argument
 */
template <typename Traits, typename KdTreeDenseBuffers, typename StaticKdTreeBuffers>
void deepCopyBuffersToDevice(sycl::queue& queue, const KdTreeDenseBuffers& hostBuffers,
                             StaticKdTreeBuffers& hostBuffersHoldingDevicePointers)
{
    using DataPoint = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using IndexType = typename Traits::IndexType; ///< Type used to index points into the PointContainer
    using NodeType  = typename Traits::NodeType;  ///< Type of nodes used inside the KdTree

    // Assign buffer sizes
    hostBuffersHoldingDevicePointers.points_size  = hostBuffers.points_size;
    hostBuffersHoldingDevicePointers.nodes_size   = hostBuffers.nodes_size;
    hostBuffersHoldingDevicePointers.indices_size = hostBuffers.indices_size;

    hostBuffersHoldingDevicePointers.points =
        sycl::malloc_device<DataPoint>(hostBuffers.points_size * sizeof(DataPoint), queue);
    hostBuffersHoldingDevicePointers.nodes =
        sycl::malloc_device<NodeType>(hostBuffers.nodes_size * sizeof(NodeType), queue);
    hostBuffersHoldingDevicePointers.indices =
        sycl::malloc_device<IndexType>(hostBuffers.indices_size * sizeof(IndexType), queue);

    queue.memcpy(hostBuffersHoldingDevicePointers.points, hostBuffers.points.data(),
                 hostBuffers.points_size * sizeof(DataPoint));
    queue.memcpy(hostBuffersHoldingDevicePointers.nodes, hostBuffers.nodes.data(),
                 hostBuffers.nodes_size * sizeof(NodeType));
    queue.memcpy(hostBuffersHoldingDevicePointers.indices, hostBuffers.indices.data(),
                 hostBuffers.indices_size * sizeof(IndexType));
    queue.wait();
}

template <typename Buffers>
void freeBuffersOnDevice(sycl::queue& queue, const Buffers& hostBuffersHoldingDevicePointers)
{
    sycl::free(hostBuffersHoldingDevicePointers.points, queue);
    sycl::free(hostBuffersHoldingDevicePointers.nodes, queue);
    sycl::free(hostBuffersHoldingDevicePointers.indices, queue);
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
