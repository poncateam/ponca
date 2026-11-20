/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <type_traits>
#include <variant>
#include <vector>

namespace nb = nanobind;

#include <Ponca/SpatialPartitioning>
#include "../Common/pypoint.h"

/**
 * \brief Abstracts away kdtree for python binding
 *
 * This function aims at providing some kind of runtime dispatch
 * between Sparse and Dense KDTree, on the CPU or on the GPU.
 *
 * Trees are always built on the CPU. toDevice copies a tree to another device, along
 * with the points of its cloud: points only hold pointers to the arrays of the cloud,
 * so the cloud is copied to the device, and the points of the tree are rebuilt from it.
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
     * \param cloud The pointcloud, which must be stored on the CPU
     * \param dense Storage type
     */
    PyKDTree(const PointCloud& cloud, bool dense = true) : m_cloud(cloud)
    {
        RequireHost(cloud.device_type(), "KDTree can only be built from arrays stored on the CPU");

        // The built tree is temporary: its buffers are copied into buffers owned by this object
        std::visit([this](const auto& tree) { StoreHostBuffers(tree.buffers()); }, Build(cloud.view(), dense));
    }

    /**
     * \brief Constructor from a PyVectorArray
     *
     * \param cloud The pointcloud
     * \param dense Storage type
     */
    PyKDTree(const PyVectorArray<Point>& cloud, bool dense = true) : PyKDTree(PointCloud(cloud), dense) {}

    /**
     * \brief Copies the tree, and the point cloud it was built from, to a device
     *
     * When the tree is already stored on this device, the returned tree shares its buffers.
     *
     * \param device Device type, which must match nanobind.
     * \param deviceId Index of the device
     */
    PyKDTree toDevice(int device, int deviceId = 0) const
    {
        if (device == device_type() && deviceId == device_id())
            return *this; // Already there, buffers are shared

        return CopyBuffersTo(m_cloud.toDevice(device, deviceId));
    }

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

    // Not exposed: only used by toDevice, which fills every member
    PyKDTree() = default;

    /**
     * \brief Copies the buffers of the tree to the device `deviceCloud` is stored on
     *
     * Points of the tree only hold pointers to the arrays of the point cloud, so they can not
     * be copied bytewise: they are rebuilt from `deviceCloud`, which must be a copy of m_cloud.
     * Nodes and indices are copied bytewise.
     *
     * \param deviceCloud Copy of m_cloud, stored on the target device
     */
    PyKDTree CopyBuffersTo(const PointCloud& deviceCloud) const
    {
        const int device   = deviceCloud.device_type();
        const int deviceId = deviceCloud.device_id();

        PyKDTree deviceTree;
        deviceTree.m_cloud = deviceCloud; // Keeps the device arrays alive, and gives the device of the tree

        // Points are built on the host from the device cloud, so the pointers they hold are valid on the device
        const auto cloudView = deviceCloud.view();
        std::vector<Point> points;
        points.reserve(cloudView.size());
        for (size_t i = 0; i < cloudView.size(); ++i)
            points.push_back(cloudView[i]);

        deviceTree.m_points  = DeviceBuffer<Point>(points, device, deviceId);
        deviceTree.m_nodes   = m_nodes.toDevice(device, deviceId);
        deviceTree.m_indices = m_indices.toDevice(device, deviceId);
        deviceTree.m_buffers = typename View::Buffers(
            deviceTree.m_points.data(), deviceTree.m_nodes.data(), deviceTree.m_indices.data(),
            deviceTree.m_points.size(), deviceTree.m_nodes.size(), deviceTree.m_indices.size()
        );
        return deviceTree;
    }

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

    /**
     * \brief Copies the buffers of a tree built on the CPU into the buffers of this object
     */
    template <typename Buffers>
    void StoreHostBuffers(const Buffers& buffers)
    {
        constexpr int cpu = nb::device::cpu::value;
        m_points          = DeviceBuffer<Point>(buffers.points, cpu);
        m_nodes           = DeviceBuffer<NodeType>(buffers.nodes, cpu);
        m_indices         = DeviceBuffer<IndexType>(buffers.indices, cpu);
        m_buffers         = typename View::Buffers(
            m_points.data(), m_nodes.data(), m_indices.data(),
            m_points.size(), m_nodes.size(), m_indices.size()
        );
    }

    // Keeps the arrays alive as the tree references them
    PointCloud m_cloud;

    // Buffers of the tree, on the device of the cloud
    DeviceBuffer<Point> m_points;
    DeviceBuffer<NodeType> m_nodes;
    DeviceBuffer<IndexType> m_indices;

    // Pointers to the buffers above. As buffers are reference counted, these remain valid in copies.
    typename View::Buffers m_buffers;
};
