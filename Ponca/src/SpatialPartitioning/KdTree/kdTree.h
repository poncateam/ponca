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

    using Vector3Array = std::vector< Eigen::Matrix<Scalar, 3, 1>>;
    using Aabb = Eigen::AlignedBox<Scalar, 3>;

    inline KdTree();
    template <typename Range>
    inline KdTree(const Range& points);
	template<class DataPoint>
    inline KdTree(std::shared_ptr<Vector3Array>& points);
	template<class DataPoint>
    inline KdTree(std::shared_ptr<Vector3Array>& points, const std::vector<int>& sampling);

    inline void clear();
    inline void build(std::shared_ptr<Vector3Array>& points);
    inline void build(std::shared_ptr<Vector3Array>& points, const std::vector<int>& sampling);
    inline void rebuild(const std::vector<int>& sampling);

    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline size_t size() const;

    inline const Vector3Array& point_data() const;
    inline       Vector3Array& point_data();

    inline const std::shared_ptr<Vector3Array>& point_ptr() const;
    inline       std::shared_ptr<Vector3Array>& point_ptr();

    inline const std::vector<KdTreeNode<Scalar>>& node_data() const;
    inline       std::vector<KdTreeNode<Scalar>>& node_data();

    inline const std::vector<int>& index_data() const;
    inline       std::vector<int>& index_data();

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
	inline KdTreeRangeIndexQuery<Scalar> range_neighbors(int index, Scalar r) const;

    // Data --------------------------------------------------------------------
protected:
    std::shared_ptr<Vector3Array>            m_points;
    std::shared_ptr<std::vector<KdTreeNode<Scalar>>> m_nodes;
    std::shared_ptr<std::vector<int>>        m_indices;

    int m_min_cell_size;
};

} // namespace Ponca

#include "./kdTree.hpp"
