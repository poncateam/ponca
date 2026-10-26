/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <optional>
#include <type_traits>
#include <variant>

namespace nb = nanobind;

#include <Ponca/SpatialPartitioning>
#include "../Common/pypoint.h"

/**
 * \brief Abstracts away kdtree for python binding
 *
 * This function aims at providing some kind of runtime dispatch
 * between Sparse and Dense KDTree, on the CPU or on the GPU.
 *
 * Trees are always built on the CPU, then their buffers are copied to the device
 * the point cloud is stored on.
 */
template <typename PointCloud>
struct PyKDTree
{
    using Point  = typename PointCloud::Point;
    using Dense  = Ponca::KdTreeDense<Point>;
    using Sparse = Ponca::KdTreeSparse<Point>;

    /**
     * \brief Kdtree type used for queries and computations, whatever the device
     *
     * Its buffers are raw pointers to the memory of the device the tree is stored on,
     * so it can be used in device code.
     */
    using View = Ponca::StaticKdTreeBase<Ponca::KdTreePointerTraits<Point>>;

    /**
     * \brief Constructor from a pointcloud
     *
     * \param cloud The pointcloud
     * \param dense Storage type
     */
    PyKDTree(const PointCloud& cloud, bool dense = true) : m_cloud(cloud)
    {
        const int device   = cloud.device_type();
        const int deviceId = cloud.device_id();

#ifdef __CUDACC__
        // Copies and allocations are performed on the current device
        // optionnal avoid constructors.
        std::optional<CudaDeviceGuard> guard;
        if (device == nb::device::cuda::value)
            guard.emplace(deviceId);
#endif

        // Only positions are required to build the tree, on the host
        std::optional<PointCloud> hostCloud;
        switch (device)
        {
        case nb::device::cpu::value:
            hostCloud = m_cloud;
            break;
#ifdef __CUDACC__
        case nb::device::cuda::value:
            hostCloud.emplace(CopyToHost(m_cloud.positions()));
            break;
#endif
        default:
            throw std::runtime_error("Unsupported device for KdTree: " + std::to_string(device));
        }

        std::visit(
            [&](const auto& tree) {
                const auto& buffers = tree.buffers();

                // Nodes and indices do not reference the point data, they can be copied as is.
                // Trees store points in the input order: the i-th point of the tree is the i-th
                // point of the cloud.
                const auto cloudView = m_cloud.view();
                std::vector<Point> points;
                points.reserve(cloudView.size());
                for (size_t i = 0; i < cloudView.size(); ++i)
                    points.push_back(cloudView[i]);

                m_points  = CopyToDevice(points, device, deviceId);
                m_nodes   = CopyToDevice(buffers.nodes, device, deviceId);
                m_indices = CopyToDevice(buffers.indices, device, deviceId);
                m_buffers = typename View::Buffers(reinterpret_cast<Point*>(m_points.data()),
                                                   reinterpret_cast<NodeType*>(m_nodes.data()),
                                                   reinterpret_cast<IndexType*>(m_indices.data()), buffers.points_size,
                                                   buffers.nodes_size, buffers.indices_size);
            },
            Build(hostCloud->view(), dense));
    }

    /**
     * \brief Constructor from a PyVectorArray
     *
     * \param cloud The pointcloud
     * \param dense Storage type
     */
    PyKDTree(const PyVectorArray<Point>& cloud, bool dense = true) : PyKDTree(PointCloud(cloud), dense) {}

    auto device_type() const { return m_cloud.device_type(); }

    auto device_id() const { return m_cloud.device_id(); }

    /**
     * \brief Returns the tree, usable on the device it is stored on
     */
    View view() const
    {
        auto buffers = m_buffers;
        return View(buffers);
    }

private:
    using NodeType  = typename View::NodeType;
    using IndexType = typename View::IndexType;

    /**
     * \brief Build a tree on the CPU
     *
     * \param cloud Points stored on the CPU
     * \param dense Storage type
     */
    static std::variant<Dense, Sparse> Build(const typename PointCloud::View& cloud, bool dense)
    {
        if (dense)
            return Dense(cloud);
        return Sparse(cloud);
    }

    // Keeps the arrays alive as the tree references them
    PointCloud m_cloud;

    // Buffers of the tree, on the device of the cloud
    nb::ndarray<uint8_t, nb::array_api> m_points;
    nb::ndarray<uint8_t, nb::array_api> m_nodes;
    nb::ndarray<uint8_t, nb::array_api> m_indices;

    // Pointers to the buffers above. As arrays are reference counted, these remain valid in copies.
    typename View::Buffers m_buffers;
};
