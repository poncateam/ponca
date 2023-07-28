/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./knnGraphTraits.h"

#include "Query/knnGraphQuery.h"
#include "Query/knnGraphRangeQuery.h"

#include "../KdTree/kdTree.h"

#include <memory>

namespace Ponca {

template <typename Traits> class KnnGraphBase;

/*!
 * \brief Public interface for KnnGraph datastructure.
 *
 * Provides default implementation of the KnnGraph
 *
 * \see KnnGraphDefaultTraits for the default trait interface documentation.
 * \see KnnGraphBase for complete API
 */
template <typename DataPoint>
using KnnGraph = KnnGraphBase<KnnGraphDefaultTraits<DataPoint>>;

/*!
 * \brief Customizable base class for KnnGraph datastructure
 *
 * \see Ponca::KnGraph
 *
 * \tparam Traits Traits type providing the types and constants used by the KnnGraph. Must have the
 * same interface as the default traits type.
 *
 * \see KnnGraphDefaultTraits for the trait interface documentation.
 *
 */
template <typename Traits> class KnnGraphBase
{
public:
    using DataPoint  = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using Scalar     = typename DataPoint::Scalar; ///< Scalar given by user via DataPoint
    using VectorType = typename DataPoint::VectorType; ///< VectorType given by user via DataPoint

    using IndexType      = typename Traits::IndexType;
    using PointContainer = typename Traits::PointContainer; ///< Container for DataPoint used inside the KdTree
    using IndexContainer = typename Traits::IndexContainer; ///< Container for indices used inside the KdTree

    using KNearestIndexQuery = KnnGraphQuery<Traits>;
    using RangeIndexQuery    = KnnGraphRangeQuery<Traits>;

    // knnGraph ----------------------------------------------------------------
public:
    /// \brief Build a KnnGraph from a KdTree
    /// \warning Stores a const reference to kdtree.point_data()
    /// \warning KdTreeTraits compatibility is checked with static assertion
    template<typename KdTreeTraits>
    inline KnnGraphBase(const KdTreeBase<KdTreeTraits>& kdtree, int k = 6)
            : m_k(k), m_points(kdtree.point_data())
    {
        static_assert( std::is_same_v<typename Traits::DataPoint, typename KdTreeTraits::DataPoint>,
                       "KdTreeTraits::DataPoint is not equal to Traits::DataPoint" );
        static_assert( std::is_same_v<typename Traits::PointContainer, typename KdTreeTraits::PointContainer>,
                       "KdTreeTraits::PointContainer is not equal to Traits::PointContainer" );

        const int size = kdtree.index_count();
        m_indices.resize(size * m_k, -1);

#pragma omp parallel for shared(kdtree, size) default(none)
        for(int i=0; i<size; ++i)
        {
            int j = 0;
            for(int n : kdtree.k_nearest_neighbors(typename KdTreeTraits::IndexType(i),
                                                   typename KdTreeTraits::IndexType(m_k)))
            {
                m_indices[i * m_k + j] = n;
                ++j;
            }
        }
    }

    // Query -------------------------------------------------------------------
public:
    inline KNearestIndexQuery k_nearest_neighbors(int index) const{
        return KnnGraphQuery(this, index);
    }

    inline RangeIndexQuery    range_neighbors(int index, Scalar r) const{
        return knnGraphRangeQuery(this, r, index);
    }

    /// \brief Get id of the i-th neighbor of a vertex
    /// \param index Vertex id
    /// \param i Neighbor id
    inline int k_neighbor(int index, int i) const{
        return m_indices->operator[](index * m_k + i);
    }

    // Accessors ---------------------------------------------------------------
public:
    inline int k() const { return m_k; }
    inline int size() const { return m_points.size(); }

    inline const PointContainer& point_data() const { return m_points; };
    inline const IndexContainer& index_data() const { return m_indices; };

    // Data --------------------------------------------------------------------
private:
    const int m_k;
    const PointContainer& m_points;
    IndexContainer m_indices; ///< \brief Stores neighborhood relations
};

} // namespace Ponca
