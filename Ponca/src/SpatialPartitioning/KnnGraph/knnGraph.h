/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./knnGraphTraits.h"

#include "Query/knnGraphQuery.h"
#include "Query/knnGraphRangeQuery.h"

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

class KdTree;
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

    using KNearestIndexQuery = KnnGraphQuery;
    using RangeIndexQuery    = KnnGraphRangeQuery;

    // knnGraph ----------------------------------------------------------------
public:
    KnnGraphBase();
    KnnGraphBase(int k);

    void clear();
    void build(const KdTree& kdtree);
    void build(const KdTree& kdtree, int k);
    void build(const KdTree& kdtree, const IndexContainer& indices);
    void build(const KdTree& kdtree, int k, const IndexContainer& indices);


    // Query -------------------------------------------------------------------
public:
    KNearestIndexQuery k_nearest_neighbors(int index) const;
    RangeIndexQuery    range_neighbors(int index, Scalar r) const;

    int k_neighbor(int index, int i) const;

    // Empty Query -------------------------------------------------------------
public:
    RangeIndexQuery range_query(Scalar r = 0) const;

    // Accessors ---------------------------------------------------------------
public:
    int k() const;
    int size() const;

    const PointContainer& point_data() const;
          PointContainer& point_data();

    const IndexContainer& index_data() const;
          IndexContainer& index_data();

    void set_verbose(bool verbose = true);

    // Data --------------------------------------------------------------------
protected:
    int m_k;
    std::shared_ptr<PointContainer>   m_points;
    std::shared_ptr<IndexContainer> m_indices;

    bool m_verbose;
};

} // namespace Ponca
