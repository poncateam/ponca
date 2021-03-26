/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./kdTreeNode.h"

#include <Eigen/Eigen>
#include <Eigen/Geometry> // aabb

#include <memory>
#include <vector>
#include <numeric>


#include "Iterator/KdTreeRangeIndexIterator.h"
#include "Iterator/KdTreeRangePointIterator.h"

#include "Query/KdTreeNearestPointQuery.h"
#include "Query/KdTreeNearestIndexQuery.h"
#include "Query/KdTreeRangeIndexQuery.h"
#include "Query/KdTreeRangePointQuery.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
template<class DataPoint>
class KdTree
{
public:
	typedef typename DataPoint::Scalar     Scalar;  // Scalar given by user
	typedef typename DataPoint::VectorType VectorType;  // VectorType given by user

	typedef typename Eigen::AlignedBox<Scalar, DataPoint::Dim> Aabb; // Intersections

    typedef typename std::vector<VectorType> PointContainer; // Container for VectorType used inside the KdTree
    typedef typename std::vector<int> IndexContainer; // Container for indices used inside the KdTree
    typedef typename std::vector<KdTreeNode<Scalar>> NodeContainer;  // Container for nodes used inside the KdTree

    inline KdTree():
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
    };

    template<typename PointUserContainer>
    inline KdTree(const PointUserContainer& points): // PointUserContainer => Given by user, transformed to PointContainer
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
        this->build(points);
    };

    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTree(const PointUserContainer& points, const IndexUserContainer& sampling): // PointUserContainer => Given by user, transformed to PointContainer
                                                                                         // IndexUserContainer => Given by user, transformed to IndexContainer
        m_points(),
        m_nodes(),
        m_indices(),
        m_min_cell_size(64)
    {
        this->build(points, sampling);
    };

    inline void clear();

    template<typename PointUserContainer>
    inline void build(const PointUserContainer& points);  // PointUserContainer => Given by user, transformed to PointContainer


    template<typename PointUserContainer, typename IndexUserContainer>
    inline void build(const PointUserContainer& points, const IndexUserContainer& sampling); // PointUserContainer => Given by user, transformed to PointContainer
                                                                                             // IndexUserContainer => Given by user, transformed to IndexContainer


    template<typename IndexUserContainer>
    inline void rebuild(const IndexUserContainer& sampling); // IndexUserContainer => Given by user, transformed to IndexContainer


    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline int node_count() const;
    inline int index_count() const;
    inline int point_count() const;

    inline PointContainer& point_data()
    {
        return *m_points.get();
    };

    inline const PointContainer& point_data() const
    {
        return m_points;
    };

    inline const PointContainer& point_ptr() const
    {
        return m_points;
    }

    inline PointContainer& point_ptr()
    {
        return m_points;
    }

    inline const NodeContainer& node_data() const
    {
        return m_nodes;
    }

    inline NodeContainer& node_data()
    {
        return *m_nodes.get();
    }

    inline const IndexContainer& index_data() const
    {
        return m_indices;
    }

    inline IndexContainer& index_data()
    {
        return *m_indices.get();
    }

    // Parameters --------------------------------------------------------------
public:
    inline int min_cell_size() const;
    inline void set_min_cell_size(int min_cell_size);

    // Internal ----------------------------------------------------------------
public:
    inline void build_rec(int node_id, int start, int end, int level);
    inline int partition(int start, int end, int dim, Scalar value);


	// Query -------------------------------------------------------------------
public :
   /* KdTreeKNearestPointQuery<DataPoint> k_nearest_neighbors(const VectorType& point, int k) const
    {
        return KdTreeKNearestPointQuery<DataPoint>(this, k, point);
    }

    KdTreeKNearestIndexQuery<DataPoint> k_nearest_neighbors(int index, int k) const
    {
        return KdTreeKNearestIndexQuery<DataPoint>(this, k, index);
    }*/

    KdTreeNearestPointQuery<DataPoint> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<DataPoint>(this, point);
    }

    KdTreeNearestIndexQuery<DataPoint> nearest_neighbor(int index) const
    {
        return KdTreeNearestIndexQuery<DataPoint>(this, index);
    }

    KdTreeRangePointQuery<DataPoint> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<DataPoint>(this, r, point);
    }

    KdTreeRangeIndexQuery<DataPoint> range_neighbors(int index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<DataPoint>(this, r, index);
    }

	// Empty Query ------------------------------------------------------------
public:
  /*  KdTreeKNearestPointQuery<DataPoint> k_nearest_point_query(int k) const
    {
        return KdTreeKNearestPointQuery<DataPoint>(this, k);
    }

    KdTreeKNearestIndexQuery<DataPoint> k_nearest_index_query(int k) const
    {
        return KdTreeKNearestIndexQuery<DataPoint>(this, k);
    }*/

    KdTreeNearestPointQuery<DataPoint> nearest_point_query() const
    {
        return KdTreeNearestPointQuery<DataPoint>(this);
    }

    KdTreeNearestIndexQuery<DataPoint> nearest_index_query() const
    {
        return KdTreeNearestIndexQuery<DataPoint>(this);
    }

    KdTreeRangePointQuery<DataPoint> range_point_query(Scalar r) const
    {
        return KdTreeRangePointQuery<DataPoint>(this, r);
    }

    KdTreeRangeIndexQuery<Scalar> range_index_query(Scalar r) const
    {
        return KdTreeRangeIndexQuery<Scalar>(this, r);
    }

    // Data --------------------------------------------------------------------
protected:
    PointContainer m_points;
    NodeContainer m_nodes;
    IndexContainer m_indices;

    int m_min_cell_size;
};

#include "./kdTree.hpp"


} // namespace Ponca
