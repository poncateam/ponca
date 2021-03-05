/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once


#include <Eigen/Eigen>
#include <Eigen/Geometry> // aabb

#include <memory>
#include <vector>
#include <numeric>

#include "./kdTreeNode.h"
#include "./Query/KdTreeRangeIndexQuery.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
template<class DataPoint>
class KdTree
{
public:

	/*! \brief Scalar type inherited from DataPoint */
	typedef typename DataPoint::Scalar     Scalar;
	/*! \brief Vector type inherited from DataPoint */
	typedef typename DataPoint::VectorType VectorType;

    using Aabb = Eigen::AlignedBox<Scalar, 3>;

    inline KdTree(): 
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
    };
    /*template <typename Range>
    inline KdTree(const Range& points);*/
    inline KdTree(std::shared_ptr<VectorType>& points):
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
        this->build(points);
    };
    inline KdTree(std::shared_ptr<VectorType>& points, const std::vector<int>& sampling):
        m_points(nullptr),
        m_nodes(nullptr),
        m_indices(nullptr),
        m_min_cell_size(64)
    {
        this->build(points, sampling);
    };

    inline void clear();
    inline void build(std::shared_ptr<VectorType>& points);
    inline void build(std::shared_ptr<VectorType>& points, const std::vector<int>& sampling);
    inline void rebuild(const std::vector<int>& sampling);

    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline size_t size() const;

    inline const VectorType& KdTree<DataPoint>::point_data() const
    {
        return *m_points.get();
    }

    inline VectorType& KdTree<DataPoint>::point_data()
    {
        return *m_points.get();
    }

    inline const std::shared_ptr<VectorType>& KdTree<DataPoint>::point_ptr() const
    {
        return m_points;
    }

    inline std::shared_ptr<VectorType>& KdTree<DataPoint>::point_ptr()
    {
        return m_points;
    }

    inline const std::vector<Ponca::KdTreeNode<Scalar>>& KdTree<DataPoint>::node_data() const
    {
        return *m_nodes.get();
    }

    inline std::vector<Ponca::KdTreeNode<Scalar>>& KdTree<DataPoint>::node_data()
    {
        return *m_nodes.get();
    }

    inline const std::vector<int>& Ponca::KdTree<DataPoint>::index_data() const
    {
        return *m_indices.get();
    }

    inline std::vector<int>& Ponca::KdTree<DataPoint>::index_data()
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
	/*KNearestPointQuery k_nearest_neighbors(const Vector3& point, int k) const;
	KNearestIndexQuery k_nearest_neighbors(int index, int k) const;
	NearestPointQuery  nearest_neighbor(const Vector3& point) const;
	NearestIndexQuery  nearest_neighbor(int index) const;
	RangePointQuery    range_neighbors(const Vector3& point, Scalar r) const;*/
	inline KdTreeRangeIndexQuery<Scalar> range_neighbors(int index, Scalar r) const;
	
	// Empty Query ------------------------------------------------------------

	/*KNearestPointQuery k_nearest_point_query(int k = 0) const;
	KNearestIndexQuery k_nearest_index_query(int k = 0) const;
	NearestPointQuery  nearest_point_query() const;
	NearestIndexQuery  nearest_index_query() const;
	RangePointQuery    range_point_query(Scalar r = 0) const;
	RangeIndexQuery    range_index_query(Scalar r = 0) const;*/

    // Data --------------------------------------------------------------------
protected:
    std::shared_ptr<VectorType>            m_points;
    std::shared_ptr<std::vector<KdTreeNode<Scalar>>> m_nodes;
    std::shared_ptr<std::vector<int>>        m_indices;

    int m_min_cell_size;
};

#include "./kdTree.hpp"

} // namespace Ponca
