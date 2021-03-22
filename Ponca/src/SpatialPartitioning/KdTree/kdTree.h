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

#include "../query.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
template<class DataPoint> class KdTreeRangeIndexIterator;
template<class Scalar> class KdTreeRangeIndexQuery;
template<class DataPoint> class KdTreeRangePointQuery;
template<class DataPoint> class KdTreeRangePointIterator;
template<class Scalar> class RangeIndexQuery;
template<class DataPoint> class KdTreeQuery;
template<class DataPoint>
class KdTree
{


public:


	/*! \brief Scalar type inherited from DataPoint */
	typedef typename DataPoint::Scalar     Scalar;
	/*! \brief Vector type inherited from DataPoint */
	typedef typename DataPoint::VectorType VectorType;
	typedef typename DataPoint::Vector Vector;

	using Aabb = Eigen::AlignedBox<Scalar, 3>;

    inline KdTree():
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
    };

    inline KdTree(std::shared_ptr<Vector>& points):
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
        this->build(points);
    };

    inline KdTree(std::shared_ptr<Vector>& points, const std::vector<int>& sampling):
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
        this->build(points, sampling);
    };

    inline void clear();
    inline void build(std::shared_ptr<Vector>& points);
    inline void build(std::shared_ptr<Vector>& points, const std::vector<int>& sampling);
    inline void rebuild(const std::vector<int>& sampling);

    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline int node_count() const;
    inline int index_count() const;
    inline int point_count() const;

    inline Vector& point_data()
    {
        return *m_points.get();
    };

    inline const Vector& point_data() const
    {
        return *m_points.get();
    };

    inline const std::shared_ptr<Vector>& point_ptr() const
    {
        return m_points;
    }

    inline std::shared_ptr<Vector>& point_ptr()
    {
        return m_points;
    }

    inline const std::vector<KdTreeNode<Scalar>>& node_data() const
    {
        return *m_nodes.get();
    }

    inline std::vector<KdTreeNode<Scalar>>& node_data()
    {
        return *m_nodes.get();
    }

    inline const std::vector<int>& index_data() const
    {
        return *m_indices.get();
    }

    inline std::vector<int>& index_data()
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
    }

    KdTreeNearestPointQuery<DataPoint> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<DataPoint>(this, point);
    }

    KdTreeNearestIndexQuery<DataPoint> nearest_neighbor(int index) const
    {
        return KdTreeNearestIndexQuery<DataPoint>(this, index);
    }*/

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
    }

    KdTreeNearestPointQuery<DataPoint> nearest_point_query() const
    {
        return KdTreeNearestPointQuery<DataPoint>(this);
    }

    KdTreeNearestIndexQuery<DataPoint> nearest_index_query() const
    {
        return KdTreeNearestIndexQuery<DataPoint>(this);
    }*/

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
    std::shared_ptr<Vector>            m_points;
    std::shared_ptr<std::vector<KdTreeNode<Scalar>>> m_nodes;
    std::shared_ptr<std::vector<int>>        m_indices;

    int m_min_cell_size;
};

#include "./kdTree.hpp"

} // namespace Ponca
