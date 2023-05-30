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

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
///
template <typename DataPoint>
struct DefaultKdTreeAdapter
{
private:
    typedef typename DataPoint::Scalar     Scalar;
    typedef typename DataPoint::VectorType VectorType;

public:
    typedef Eigen::AlignedBox<Scalar, DataPoint::Dim> AabbType;

    typedef int DimType;
    typedef int DepthType;

    // Containers
    typedef std::vector<DataPoint> PointContainer;
    typedef std::vector<int>       IndexContainer;

    typedef std::vector<DefaultKdTreeNode<DataPoint, AabbType>> NodeContainer;

    static DimType max_dim(const VectorType& vec)
    {
        DimType dim;
        vec.maxCoeff(dim);
        return dim;
    }

    static Scalar vec_component(const VectorType& vec, DimType dim)
    {
        return vec(dim);
    }
};

///
/// \tparam DataPoint
/// \tparam Adapter
///
/// \todo Better handle sampling: do not store non-selected points (requires to store original indices
template<class _DataPoint, class Adapter = DefaultKdTreeAdapter<_DataPoint>>
class KdTree
{
public:
    typedef          _DataPoint            DataPoint;
    typedef typename DataPoint::Scalar     Scalar; // Scalar given by user
    typedef typename DataPoint::VectorType VectorType; // VectorType given by user

    typedef typename Adapter::AabbType AabbType;

    typedef typename Adapter::DimType   DimType;
    typedef typename Adapter::DepthType DepthType;

    typedef typename Adapter::PointContainer PointContainer; // Container for VectorType used inside the KdTree
    typedef typename Adapter::IndexContainer IndexContainer; // Container for indices used inside the KdTree
    typedef typename Adapter::NodeContainer  NodeContainer; // Container for nodes used inside the KdTree

    typedef typename IndexContainer::value_type IndexType;
    typedef typename NodeContainer::value_type  NodeType;
    typedef typename PointContainer::size_type  PointCountType;
    typedef typename IndexContainer::size_type  IndexCountType;
    typedef typename NodeContainer::size_type   NodeCountType;

    typedef typename NodeType::LeafSizeType LeafSizeType;

    static_assert(std::is_same<typename PointContainer::value_type, DataPoint>::value, "Point container must contain DataPoints");

    inline KdTree():
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
    };

    template<typename PointUserContainer>
    inline KdTree(PointUserContainer&& points): // PointUserContainer => Given by user, transformed to PointContainer
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        this->build(std::forward<PointUserContainer>(points));
    };

    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTree(PointUserContainer&& points, IndexUserContainer sampling): // PointUserContainer => Given by user, transformed to PointContainer
                                                                             // IndexUserContainer => Given by user, transformed to IndexContainer
        m_points(),
        m_nodes(),
        m_indices(),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling));
    };

    inline void clear();

    struct DefaultConverter{
        template <typename Input>
        inline void operator()( Input&& i, PointContainer & o ) {
            typedef typename std::remove_cv<typename std::remove_reference<Input>::type>::type InputContainer;
            if constexpr ( std::is_same<InputContainer, PointContainer>::value )
                o = i;
            else
                std::transform(i.cbegin(), i.cend(), std::back_inserter(o),
                               [](const typename InputContainer::value_type &p) -> DataPoint { return DataPoint(p); });
        }
    };

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \param points Input points
    template<typename PointUserContainer>
    inline void build(PointUserContainer&& points)
    {
        build(std::forward<PointUserContainer>(points), DefaultConverter());
    }

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling container, transformed to IndexContainer
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename Converter>
    inline void build(PointUserContainer&& points, Converter c);

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline void buildWithSampling(PointUserContainer&& points,
                                  IndexUserContainer sampling)
    {
        buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling), DefaultConverter());
    }

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam Converter
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
    inline void buildWithSampling(PointUserContainer&& points,
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

    inline void set_min_cell_size(LeafSizeType min_cell_size)
    {
        PONCA_DEBUG_ASSERT(min_cell_size > 0);
        m_min_cell_size = min_cell_size;
    }

    // Internal ----------------------------------------------------------------
public:
    inline void build_rec(NodeCountType node_id, IndexCountType start, IndexCountType end, DepthType level);
    inline IndexCountType partition(IndexCountType start, IndexCountType end, DimType dim, Scalar value);

    // Query -------------------------------------------------------------------
public :
    KdTreeKNearestPointQuery<DataPoint, Adapter> k_nearest_neighbors(const VectorType& point, int k) const
    {
        return KdTreeKNearestPointQuery<DataPoint, Adapter>(this, k, point);
    }

    KdTreeKNearestIndexQuery<DataPoint, Adapter> k_nearest_neighbors(int index, int k) const
    {
        return KdTreeKNearestIndexQuery<DataPoint, Adapter>(this, k, index);
    }

    KdTreeNearestPointQuery<DataPoint, Adapter> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<DataPoint, Adapter>(this, point);
    }

    KdTreeNearestIndexQuery<DataPoint, Adapter> nearest_neighbor(int index) const
    {
        return KdTreeNearestIndexQuery<DataPoint, Adapter>(this, index);
    }

    KdTreeRangePointQuery<DataPoint, Adapter> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<DataPoint, Adapter>(this, r, point);
    }

    KdTreeRangeIndexQuery<DataPoint, Adapter> range_neighbors(int index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<DataPoint, Adapter>(this, r, index);
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
