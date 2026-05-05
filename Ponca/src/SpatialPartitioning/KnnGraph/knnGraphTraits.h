/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <Eigen/Geometry>

namespace Ponca
{

    /*!
     * \brief The default traits type used by the kd-tree.
     */
    template <typename _DataPoint>
    struct KnnGraphDefaultTraits
    {
        enum
        {
            MAX_RANGE_NEIGHBORS_SIZE = 128 //!< The maximum number of neighbors in a range neighbors query
        };
        /*!
         * \brief The type used to store point data.
         *
         * Must provide `Scalar` and `VectorType` aliases.
         *
         * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
         * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
         * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
         */
        using DataPoint = _DataPoint;

    private:
        using Scalar     = typename DataPoint::Scalar;
        using VectorType = typename DataPoint::VectorType;

    public:
        /*!
         * \brief The type used to calculate node bounding boxes.
         *
         * Must provide `min()`, `max()`, and `center()` functions, all returning a `VectorType`.
         */
        using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

        // Containers
        using IndexType = int;
        /// \brief Type used to store the external Point container in the KnnGraph::Buffer
        using PointContainer = const std::vector<DataPoint>&;
        /// \brief Type used to store the index container in the KnnGraph::Buffer
        using IndexContainer = std::vector<IndexType>;
        /// \brief Type to be used to send the index container as function parameter
        using IndexContainerRef = IndexContainer&;

        /// \brief Provides access to the raw pointer where indices are stored
        static IndexType* getIndexRawPtr(IndexContainer& idx) { return idx.data(); }
        /// \brief Provides access to the raw pointer where indices are stored
        static const IndexType* getIndexRawPtr(const IndexContainer& idx) { return idx.data(); }
    };
    /*!
     * \brief Variant to the KnnGraph Traits type that uses pointers as internal storage instead of an STL-like
     * container.
     */
    template <typename _DataPoint>
    struct KnnGraphPointerTraits
    {
        enum
        {
            MAX_RANGE_NEIGHBORS_SIZE = 128 //!< The maximum number of neighbors in a range neighbors query
        };
        /*!
         * \brief The type used to store point data.
         *
         * Must provide `Scalar` and `VectorType` aliases.
         *
         * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
         * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
         * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
         */
        using DataPoint = _DataPoint;

    private:
        using Scalar     = typename DataPoint::Scalar;
        using VectorType = typename DataPoint::VectorType;

    public:
        /*!
         * \brief The type used to calculate node bounding boxes.
         *
         * Must provide `min()`, `max()`, and `center()` functions, all returning a `VectorType`.
         */
        using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

        // Containers
        using IndexType = int;
        /// \brief Type used to store the external Point container in the KnnGraph::Buffer
        /// Non-const to allow KnnGraph::Buffers copy and writing to other devices
        using PointContainer = DataPoint*;
        /// \brief Type used to store the index container in the KnnGraph::Buffer
        using IndexContainer = IndexType*;
        /// \brief Type to be used to send the index container as function parameter
        using IndexContainerRef = IndexContainer;

        /// \brief Provides access to the raw pointer where indices are stored
        static IndexType* getIndexRawPtr(IndexContainer& idx) { return idx; }
        /// \brief Provides access to the raw pointer where indices are stored
        static const IndexType* getIndexRawPtr(const IndexContainer& idx) { return idx; }
    };
} // namespace Ponca
