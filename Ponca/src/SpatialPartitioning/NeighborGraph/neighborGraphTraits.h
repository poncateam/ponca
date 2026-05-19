/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <Eigen/Geometry>
#include <stack>
#include <set>

#include "../../Common/Containers/hashset.h"
#include "../../Common/Containers/stack.h"

namespace Ponca
{

    /*!
     * \brief The default traits type used by the neighbor graph.
     */
    template <typename _DataPoint>
    struct NeighborGraphDefaultTraits
    {
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
        /// \brief Type used to store the external Point container in AbstractNeighborGraph::Buffers
        using PointContainer = const std::vector<DataPoint>&;
        /// \brief Type used to store the index container in the AbstractNeighborGraph::Buffers
        using IndexContainer = std::vector<IndexType>;
        /// \brief Type to be used to send the index container as function parameter
        using IndexContainerRef = IndexContainer&;

        /*! \brief A Set dynamic in memory, used by NeighborGraphRangeQuery
         *  \warning Not compatible with CUDA
         */
        using KnnGraphRangeSet = std::set<int>;
        /*! \brief A Stack dynamic in memory, used by NeighborGraphRangeQuery
         *  \warning Not compatible with CUDA
         */
        using KnnGraphRangeStack = std::stack<int>;

        /// \brief Provides access to the raw pointer where indices are stored
        PONCA_MULTIARCH static IndexType* getIndexRawPtr(IndexContainer& idx) { return idx.data(); }
        /// \brief Provides access to the raw pointer where indices are stored
        PONCA_MULTIARCH static const IndexType* getIndexRawPtr(const IndexContainer& idx) { return idx.data(); }
    };
    /*!
     * \brief Variant to the NeighborGraphDefaultTraits that uses pointers as internal storage instead of an STL-like
     * container.
     */
    template <typename _DataPoint>
    struct NeighborGraphPointerTraits
    {
        enum
        {
            K_MAX_NN = 1000 //!< The maximum number of neighbors that will be visited in a range neighbors query
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
        /// \brief Type used to store the external Point container in AbstractNeighborGraph::Buffers
        /// Non-const to allow AbstractNeighborGraph::Buffers copy and writing to other devices
        using PointContainer = DataPoint*;
        /// \brief Type used to store the index container in the AbstractNeighborGraph::Buffers
        using IndexContainer = IndexType*;
        /// \brief Type to be used to send the index container as function parameter
        using IndexContainerRef = IndexContainer;

        //! \brief A static Set used by NeighborGraphRangeQuery
        using KnnGraphRangeSet = HashSet<K_MAX_NN>;
        //! \brief A static Stack used by NeighborGraphRangeQuery
        using KnnGraphRangeStack = Stack<int, K_MAX_NN>;

        /// \brief Provides access to the raw pointer where indices are stored
        PONCA_MULTIARCH static IndexType* getIndexRawPtr(IndexContainer& idx) { return idx; }
        /// \brief Provides access to the raw pointer where indices are stored
        PONCA_MULTIARCH static const IndexType* getIndexRawPtr(const IndexContainer& idx) { return idx; }
    };
} // namespace Ponca
