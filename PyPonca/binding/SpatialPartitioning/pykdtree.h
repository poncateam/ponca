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

namespace nb = nanobind;

#include <Ponca/SpatialPartitioning>
#include "../Common/pypoint.h"

/**
 * \brief Abstracts away kdtree for python binding
 *
 * This function aims at providing some kind of runtime dispatch
 * between Sparse and Dense KDTree. In the future, it may also
 * be used for CPU/GPU trees.
 */
template <typename PointCloud>
struct PyKDTree
{
    using Point  = typename PointCloud::Point;
    using Dense  = Ponca::KdTreeDense<Point>;
    using Sparse = Ponca::KdTreeSparse<Point>;

    PyKDTree() : isDense(true) {}

    /**
     * \brief Constructor from a pointcloud
     *
     * \param cloud The pointcloud
     * \param dense Storage type
     */
    PyKDTree(const PointCloud& cloud, bool dense = true) : isDense(dense)
    {
        if (isDense)
            data = Dense(cloud);
        else
            data = Sparse(cloud);
    }

    /**
     * \brief Constructor from a PyVectorArray
     *
     * \param cloud The pointcloud
     * \param dense Storage type
     */
    PyKDTree(const PyVectorArray<Point>& cloud, bool dense = true) : isDense(dense)
    {
        PointCloud tmp(cloud);
        if (isDense)
            data = Dense(tmp);
        else
            data = Sparse(tmp);
    }

    /**
     * \brief Run a function with a dispatch depending on kdtree type
     *
     * \tparam Func Function type
     *
     * \param f The function
     */
    template <typename Func>
    void Run(Func&& f) const
    {
        if (isDense)
            f(std::get<DenseIndex>(data));
        else
            f(std::get<SparseIndex>(data));
    }

    auto device_type() const
    {
        return nb::device::cpu::value;
    }

    // Should be const but copy operator makes it impossible...
    bool isDense;
    static constexpr unsigned int DenseIndex  = 0;
    static constexpr unsigned int SparseIndex = 1;
    std::variant<Dense, Sparse> data;
};

