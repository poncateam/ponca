/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./knnGraphTraits.h"

#include "Query/knnGraphKNearestQuery.h"
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

    using KNearestIndexQuery = KnnGraphKNearestQuery<Traits>;
    using RangeIndexQuery    = KnnGraphRangeQuery<Traits>;

    friend class KnnGraphKNearestQuery<Traits>; // This type must be equal to KnnGraphBase::KNearestIndexQuery
    friend class KnnGraphRangeQuery<Traits>;    // This type must be equal to KnnGraphBase::RangeIndexQuery

    // knnGraph ----------------------------------------------------------------
public:
    /// \brief Build a KnnGraph from a KdTreeDense
    ///
    /// \param k Number of requested neighbors. Might be reduced if k is larger than the kdtree size - 1
    ///          (query point is not included in query output, thus -1)
    ///
    /// \warning Stores a const reference to kdtree.point_data()
    /// \warning KdTreeTraits compatibility is checked with static assertion
    template<typename KdTreeTraits>
    inline KnnGraphBase(const KdTreeBase<KdTreeTraits>& kdtree, int k = 6)
            : m_k(std::min(k,kdtree.sample_count()-1)),
              m_kdTreePoints(kdtree.points())
    {
        static_assert( std::is_same<typename Traits::DataPoint, typename KdTreeTraits::DataPoint>::value,
                       "KdTreeTraits::DataPoint is not equal to Traits::DataPoint" );
        static_assert( std::is_same<typename Traits::PointContainer, typename KdTreeTraits::PointContainer>::value,
                       "KdTreeTraits::PointContainer is not equal to Traits::PointContainer" );
        static_assert( std::is_same<typename Traits::IndexContainer, typename KdTreeTraits::IndexContainer>::value,
                       "KdTreeTraits::IndexContainer is not equal to Traits::IndexContainer" );

        // We need to account for the entire point set, irrespectively of the sampling. This is because the kdtree
        // (k_nearest_neighbors) return ids of the entire point set, not it sub-sampled list of ids.
        // \fixme Update API to properly handle kdtree subsampling
        const int cloudSize   = kdtree.point_count();
        {
            const int samplesSize = kdtree.sample_count();
            eigen_assert(cloudSize == samplesSize);
        }

        m_indices.resize(cloudSize * m_k, -1);

#pragma omp parallel for shared(kdtree, cloudSize) default(none)
        for(int i=0; i<cloudSize; ++i)
        {
            int j = 0;
            for(auto n : kdtree.k_nearest_neighbors(typename KdTreeTraits::IndexType(i),
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
        return KNearestIndexQuery(this, index);
    }

    inline RangeIndexQuery    range_neighbors(int index, Scalar r) const{
        return RangeIndexQuery(this, r, index);
    }

    // Accessors ---------------------------------------------------------------
public:
    /// \brief Number of neighbor per vertex
    inline int k() const { return m_k; }
    /// \brief Number of vertices in the neighborhood graph
    inline int size() const { return static_cast<int>(m_indices.size())/m_k; }

    // Data --------------------------------------------------------------------
private:
    const int m_k;
    IndexContainer m_indices; ///< \brief Stores neighborhood relations

protected: // for friends relations
    const PointContainer& m_kdTreePoints;
    inline const IndexContainer& index_data() const { return m_indices; };
};

} // namespace Ponca
