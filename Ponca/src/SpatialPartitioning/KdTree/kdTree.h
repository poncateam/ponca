/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./kdTreeTraits.h"

#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Geometry> // aabb

#include "../../Common/Assert.h"

#include "Query/kdTreeNearestPointQuery.h"
#include "Query/kdTreeNearestIndexQuery.h"
#include "Query/kdTreeKNearestPointQuery.h"
#include "Query/kdTreeKNearestIndexQuery.h"
#include "Query/kdTreeRangeIndexQuery.h"
#include "Query/kdTreeRangePointQuery.h"

namespace Ponca {
template <typename Traits> class KdTreeBase;

///
/// \tparam DataPoint
template <typename DataPoint>
using KdTree = KdTreeBase<KdTreeDefaultTraits<DataPoint>>;

///
/// \tparam Traits
///
/// \todo Finish documentation
/// The specified traits type must have the following interface:
///
/// DataPoint
///
/// AabbType
///
/// IndexType
/// PointContainer
/// IndexContainer
/// NodeContainer
///
/// MAX_DEPTH
///
/// PointContainer::value_type must match DataPoint
/// IndexContainer::value_type must match IndexType
///
/// \todo Better handle sampling: do not store non-selected points (requires to store original indices)
template <typename Traits>
class KdTreeBase
{
public:
    using DataPoint  = typename Traits::DataPoint;
    using Scalar     = typename DataPoint::Scalar; // Scalar given by user
    using VectorType = typename DataPoint::VectorType; // VectorType given by user
    using AabbType   = typename Traits::AabbType;

    using IndexType      = typename Traits::IndexType;
    using PointContainer = typename Traits::PointContainer; // Container for DataPoint used inside the KdTree
    using IndexContainer = typename Traits::IndexContainer; // Container for indices used inside the KdTree
    using NodeContainer  = typename Traits::NodeContainer; // Container for nodes used inside the KdTree

    using NodeType      = typename NodeContainer::value_type;
    using NodeCountType = typename NodeContainer::size_type;
    using LeafSizeType  = typename NodeType::LeafSizeType;

    static_assert(std::is_same<typename PointContainer::value_type, DataPoint>::value,
        "PointContainer must contain DataPoints");
    
    // Queries use a value of -1 for invalid indices
    static_assert(std::is_signed<IndexType>::value, "Index type must be signed");
    static_assert(std::is_same<typename IndexContainer::value_type, IndexType>::value, "Index type mismatch");

    static_assert(Traits::MAX_DEPTH > 0, "Max depth must be strictly positive");

    inline KdTreeBase():
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
    };

    template<typename PointUserContainer>
    inline KdTreeBase(PointUserContainer points): // PointUserContainer => Given by user, transformed to PointContainer
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        this->build(std::move(points));
    };

    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTreeBase(PointUserContainer points, IndexUserContainer sampling): // PointUserContainer => Given by user, transformed to PointContainer
                                                                               // IndexUserContainer => Given by user, transformed to IndexContainer
        m_points(),
        m_nodes(),
        m_indices(),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        buildWithSampling(std::move(points), std::move(sampling));
    };

    inline void clear();

    struct DefaultConverter
    {
        template <typename Input>
        inline void operator()(Input i, PointContainer& o)
        {
            if constexpr (std::is_same<Input, PointContainer>::value)
                o = std::move(i);
            else
                std::transform(i.cbegin(), i.cend(), std::back_inserter(o),
                               [](const typename Input::value_type &p) -> DataPoint { return DataPoint(p); });
        }
    };

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \param points Input points
    template<typename PointUserContainer>
    inline void build(PointUserContainer points)
    {
        build(std::move(points), DefaultConverter());
    }

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling container, transformed to IndexContainer
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename Converter>
    inline void build(PointUserContainer points, Converter c);

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline void buildWithSampling(PointUserContainer points,
                                  IndexUserContainer sampling)
    {
        buildWithSampling(std::move(points), std::move(sampling), DefaultConverter());
    }

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam Converter
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
    inline void buildWithSampling(PointUserContainer points,
                                  IndexUserContainer sampling,
                                  Converter c);


    template<typename IndexUserContainer>
    inline void rebuild(IndexUserContainer sampling); // IndexUserContainer => Given by user, transformed to IndexContainer


    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline NodeCountType node_count() const
    {
        return m_nodes.size();
    }

    inline IndexType index_count() const
    {
        return (IndexType)m_indices.size();
    }

    inline IndexType point_count() const
    {
        return (IndexType)m_points.size();
    }

    inline NodeCountType leaf_count() const
    {
        return m_leaf_count;
    }

    inline PointContainer& point_data()
    {
        return m_points;
    };

    inline const PointContainer& point_data() const
    {
        return m_points;
    };

    inline const NodeContainer& node_data() const
    {
        return m_nodes;
    }

    inline NodeContainer& node_data()
    {
        return m_nodes;
    }

    inline const IndexContainer& index_data() const
    {
        return m_indices;
    }

    inline IndexContainer& index_data()
    {
        return m_indices;
    }

    // Parameters --------------------------------------------------------------
public:
    inline LeafSizeType min_cell_size() const
    {
        return m_min_cell_size;
    }

    inline void set_min_cell_size(LeafSizeType min_cell_size)
    {
        PONCA_DEBUG_ASSERT(min_cell_size > 0);
        m_min_cell_size = min_cell_size;
    }

    // Internal ----------------------------------------------------------------
public:
    inline void build_rec(NodeCountType node_id, IndexType start, IndexType end, int level);
    inline IndexType partition(IndexType start, IndexType end, int dim, Scalar value);

    // Query -------------------------------------------------------------------
public :
    KdTreeKNearestPointQuery<Traits> k_nearest_neighbors(const VectorType& point, IndexType k) const
    {
        return KdTreeKNearestPointQuery<Traits>(this, k, point);
    }

    KdTreeKNearestIndexQuery<Traits> k_nearest_neighbors(IndexType index, IndexType k) const
    {
        return KdTreeKNearestIndexQuery<Traits>(this, k, index);
    }

    KdTreeNearestPointQuery<Traits> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<Traits>(this, point);
    }

    KdTreeNearestIndexQuery<Traits> nearest_neighbor(IndexType index) const
    {
        return KdTreeNearestIndexQuery<Traits>(this, index);
    }

    KdTreeRangePointQuery<Traits> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<Traits>(this, r, point);
    }

    KdTreeRangeIndexQuery<Traits> range_neighbors(IndexType index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<Traits>(this, r, index);
    }
    

    // Data --------------------------------------------------------------------
protected:
    PointContainer m_points;
    NodeContainer m_nodes;
    IndexContainer m_indices;

    LeafSizeType m_min_cell_size;
    NodeCountType m_leaf_count; ///< Number of leaves in the Kdtree (computed during construction)
};

#include "./kdTree.hpp"
} // namespace Ponca
