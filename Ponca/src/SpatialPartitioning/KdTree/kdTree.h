/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./kdTreeNode.h"

#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include "../../Common/Assert.h"

#include "Query/kdTreeNearestPointQuery.h"
#include "Query/kdTreeNearestIndexQuery.h"
#include "Query/kdTreeKNearestPointQuery.h"
#include "Query/kdTreeKNearestIndexQuery.h"
#include "Query/kdTreeRangeIndexQuery.h"
#include "Query/kdTreeRangePointQuery.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
///
template <typename DataPoint>
struct DefaultKdTreeCompatibility
{
    typedef typename DataPoint::Scalar     Scalar;
    typedef typename DataPoint::VectorType VectorType;

    typedef int DimType;
    typedef int DepthType;

    // Containers
    typedef std::vector<DataPoint> PointContainer;
    typedef std::vector<int>       IndexContainer;

    typedef std::vector<DefaultKdTreeNode<DataPoint>> NodeContainer;
};

///
/// \tparam DataPoint
///
/// \todo Better handle sampling: do not store non-selected points (requires to store original indices
template<class _DataPoint,
         class Compatibility = DefaultKdTreeCompatibility<_DataPoint>>
class KdTree
{
public:
    typedef          _DataPoint                DataPoint;
    typedef typename Compatibility::Scalar     Scalar; // Scalar given by user
    typedef typename Compatibility::VectorType VectorType; // VectorType given by user

    typedef typename Compatibility::DimType   DimType;
    typedef typename Compatibility::DepthType DepthType;

    typedef typename Compatibility::PointContainer PointContainer; // Container for VectorType used inside the KdTree
    typedef typename Compatibility::IndexContainer IndexContainer; // Container for indices used inside the KdTree
    typedef typename Compatibility::NodeContainer  NodeContainer; // Container for nodes used inside the KdTree

    typedef typename IndexContainer::value_type IndexType;
    typedef typename NodeContainer::value_type  NodeType;
    typedef typename PointContainer::size_type  PointCountType;
    typedef typename IndexContainer::size_type  IndexCountType;
    typedef typename NodeContainer::size_type   NodeCountType;

    typedef typename NodeType::LeafSizeType LeafSizeType;

    static_assert(std::is_same_v<typename PointContainer::value_type, DataPoint>, "Point container must contain DataPoints");

    inline KdTree():
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
    };

    template<typename PointUserContainer>
    inline KdTree(const PointUserContainer& points): // PointUserContainer => Given by user, transformed to PointContainer
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        this->build(points);
    };

    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTree(const PointUserContainer& points, const IndexUserContainer& sampling): // PointUserContainer => Given by user, transformed to PointContainer
                                                                                         // IndexUserContainer => Given by user, transformed to IndexContainer
        m_points(),
        m_nodes(),
        m_indices(),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        buildWithSampling(points, sampling);
    };

    inline void clear();

    struct DefaultConverter{
        template <typename Input>
        inline void operator()( const Input& i, PointContainer & o ) {
            if constexpr ( std::is_same<Input, PointContainer>::value )
                o = i;
            else
                std::transform(i.cbegin(), i.cend(), std::back_inserter(o),
                               [](const typename Input::value_type &p) -> DataPoint { return DataPoint(p); });
        }
    };

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \param points Input points
    template<typename PointUserContainer>
    inline void build(const PointUserContainer& points) { build(points, DefaultConverter()); }
    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling container, transformed to IndexContainer
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename Converter>
    inline void build(const PointUserContainer& points, Converter c);

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline void buildWithSampling(const PointUserContainer& points,
                                  const IndexUserContainer& sampling)
                                  { buildWithSampling(points, sampling, DefaultConverter());}

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam Converter
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
    inline void buildWithSampling(const PointUserContainer& points,
                                  const IndexUserContainer& sampling,
                                  Converter c);


    template<typename IndexUserContainer>
    inline void rebuild(const IndexUserContainer& sampling); // IndexUserContainer => Given by user, transformed to IndexContainer


    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline NodeCountType node_count() const
    {
        return m_nodes.size();
    }

    inline IndexCountType index_count() const
    {
        return m_indices.size();
    }

    inline PointCountType point_count() const
    {
        return m_points.size();
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

    inline void set_min_cell_size(int min_cell_size)
    {
        m_min_cell_size = min_cell_size;
    }

    // Internal ----------------------------------------------------------------
public:
    inline void build_rec(NodeCountType node_id, IndexCountType start, IndexCountType end, DepthType level);
    inline IndexCountType partition(IndexCountType start, IndexCountType end, DimType dim, Scalar value);


    // Query -------------------------------------------------------------------
public :
    KdTreeKNearestPointQuery<DataPoint, Compatibility> k_nearest_neighbors(const VectorType& point, int k) const
    {
        return KdTreeKNearestPointQuery<DataPoint, Compatibility>(this, k, point);
    }

    KdTreeKNearestIndexQuery<DataPoint, Compatibility> k_nearest_neighbors(int index, int k) const
    {
        return KdTreeKNearestIndexQuery<DataPoint, Compatibility>(this, k, index);
    }

    KdTreeNearestPointQuery<DataPoint, Compatibility> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<DataPoint, Compatibility>(this, point);
    }

    KdTreeNearestIndexQuery<DataPoint, Compatibility> nearest_neighbor(int index) const
    {
        return KdTreeNearestIndexQuery<DataPoint, Compatibility>(this, index);
    }

    KdTreeRangePointQuery<DataPoint, Compatibility> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<DataPoint, Compatibility>(this, r, point);
    }

    KdTreeRangeIndexQuery<DataPoint, Compatibility> range_neighbors(int index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<DataPoint, Compatibility>(this, r, index);
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
