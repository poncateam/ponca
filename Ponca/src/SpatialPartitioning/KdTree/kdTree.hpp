#include "kdTree.h"
/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/



#include <numeric>

// KdTree ----------------------------------------------------------------------


template<class DataPoint>
inline Ponca::KdTree<DataPoint>::KdTree() : 
	m_points(nullptr),
	m_nodes(nullptr),
	m_indices(nullptr),
	m_min_cell_size(64)
{
}

template<class DataPoint>
inline Ponca::KdTree<DataPoint>::KdTree(std::shared_ptr<DataPoint::VectorType>& points) :
	m_points(nullptr),
	m_nodes(nullptr),
	m_indices(nullptr),
	m_min_cell_size(64)
{
	this->build(points);
}

template<class DataPoint>
inline Ponca::KdTree<DataPoint>::KdTree(std::shared_ptr<DataPoint::VectorType>& points, const std::vector<int>& sampling) :
	m_points(nullptr),
	m_nodes(nullptr),
	m_indices(nullptr),
	m_min_cell_size(64)
{
	this->build(points, sampling);
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::clear()
{
	m_points = nullptr;
	m_nodes   = nullptr;
	m_indices = nullptr;
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::build(std::shared_ptr<DataPoint::VectorType>& points)
{
	std::vector<int> ids;
	std::iota(ids.begin(), ids.end(), 0);
	this->build( points, ids);
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::build(std::shared_ptr<DataPoint::VectorType>& points, const std::vector<int>& sampling)
{
	this->clear();
	
	m_points = points;
	
	m_nodes = std::make_shared<std::vector<KdTreeNode>>();
	m_nodes->reserve(4 * m_points->size() / m_min_cell_size);
	m_nodes->emplace_back();
	m_nodes->back().leaf = false;
	
	m_indices = std::make_shared<std::vector<int>>(sampling);
	
	int end = static_cast<int>(m_indices->size());
	
	this->build_rec(0, 0, end, 1);
	
//    PCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::rebuild(const std::vector<int>& sampling)
{
	//    PCA_DEBUG_ASSERT(sampling.size() <= m_points->size());

    m_nodes->clear();
    m_nodes->emplace_back();
    m_nodes->back().leaf = false;

    *m_indices = sampling;

    int end = static_cast<int>(m_indices->size());
    this->build_rec(0, 0, end, 1);

//    PCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
inline bool Ponca::KdTree<DataPoint>::valid() const
{
	if (m_points == nullptr)
		return m_nodes == nullptr && m_indices == nullptr;
		
	if(m_nodes == nullptr || m_indices == nullptr)
	{
		//PCA_DEBUG_ERROR;
		return false;
	}
		
	if(m_points->size() < m_indices->size())
	{
		//PCA_DEBUG_ERROR;
		return false;
	}
		
	std::vector<bool> b(m_points->size(), false);
	for(int idx : *m_indices.get())
	{
		if(idx < 0 || int(m_points->size()) <= idx || b[idx])
		{
		    //PCA_DEBUG_ERROR;
		    return false;
		}
		b[idx] = true;
	}
		
	for(size_t n=0; n<m_nodes->size(); ++n)
	{
		const KdTreeNode& node = m_nodes->operator [](n);
		if(node.leaf)
		{
		    if(m_indices->size() <= node.start || m_indices->size() < node.start+node.size)
		    {
		        //PCA_DEBUG_ERROR;
		        return false;
		    }
		}
		else
		{
		    if(node.dim < 0 || 2 < node.dim)
		    {
		        //PCA_DEBUG_ERROR;
		        return false;
		    }
		    if(m_nodes->size() <= node.firstChildId || m_nodes->size() <= node.firstChildId+1u)
		    {
		        //PCA_DEBUG_ERROR;
		        return false;
		    }
		}
	}
		
	return true;
}

template<class DataPoint>
inline std::string Ponca::KdTree<DataPoint>::to_string() const
{
	if (!m_indices) return "";
	
	std::stringstream str;
	str << "indices (" << m_indices->size() << ") :\n";
	for(size_t i=0; i<m_indices->size(); ++i)
	{
	    str << "  " << i << ": " << m_indices->operator[](i) << "\n";
	}
	str << "nodes (" << m_nodes->size() << ") :\n";
	for(size_t n=0; n<m_nodes->size(); ++n)
	{
	    const KdTreeNode& node = m_nodes->operator[](n);
	    if(node.leaf)
	    {
	        int end = node.start + node.size;
	        str << "  leaf: start=" << node.start << " end=" << end << " (size=" << node.size << ")\n";
	    }
	    else
	    {
	        str << "  node: dim=" << node.dim << " split=" << node.splitValue << " child=" << node.firstChildId << "\n";
	    }
	}
	return str.str();
}

template<class DataPoint>
inline size_t Ponca::KdTree<DataPoint>::size() const
{
	return size_t();
}


template<class DataPoint>
inline const DataPoint::VectorType & Ponca::KdTree<DataPoint>::point_data() const
{
	return *m_points.get();
}

template<class DataPoint>
inline  DataPoint::VectorType & Ponca::KdTree<DataPoint>::point_data()
{
	return *m_points.get();
}

template<class DataPoint>
inline const std::shared_ptr< DataPoint::VectorType>& Ponca::KdTree<DataPoint>::point_ptr() const
{
	return m_points;
}

template<class DataPoint>
inline std::shared_ptr<DataPoint::VectorType>& Ponca::KdTree<DataPoint>::point_ptr()
{
	return m_points;
}

template<class DataPoint>
inline const std::vector<Ponca::KdTreeNode<DataPoint::Scalar>>& Ponca::KdTree<DataPoint>::node_data() const
{
	return *m_nodes.get();
}

template<class DataPoint>
inline std::vector<Ponca::KdTreeNode<DataPoint::Scalar>>& Ponca::KdTree<DataPoint>::node_data()
{
	return *m_nodes.get();
}

template<class DataPoint>
inline const std::vector<int>& Ponca::KdTree<DataPoint>::index_data() const
{
	return *m_indices.get();
}

template<class DataPoint>
inline std::vector<int>& Ponca::KdTree<DataPoint>::index_data()
{
	return *m_indices.get();
}

template<class DataPoint>
inline int Ponca::KdTree<DataPoint>::min_cell_size() const
{
	return m_min_cell_size;
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::set_min_cell_size(int min_cell_size)
{
	m_min_cell_size = min_cell_size;
}

template<class DataPoint>
inline void Ponca::KdTree<DataPoint>::build_rec(int node_id, int start, int end, int level)
{
	auto& nodes = *m_nodes.get();
	const auto& points  = *m_points.get();
	const auto& indices = *m_indices.get();
	
	KdTreeNode& node = nodes[node_id];
	Aabb aabb;
	for(int i=start; i<end; ++i)
	    aabb.extend(points[indices[i]]);
	
	int dim;
	(Scalar(0.5)*(aabb.max()-aabb.min())).maxCoeff(&dim);
	
	node.dim = dim;
	node.splitValue = aabb.center()(dim);
	
	int midId = this->partition(start, end, dim, node.splitValue);
	node.firstChildId = nodes.size();
	
	{
	    KdTreeNode n;
	    n.size = 0;
	    nodes.push_back(n);
	    nodes.push_back(n);
	}
	{
	    // left child
	    int childId = nodes[node_id].firstChildId;
	    KdTreeNode& child = nodes[childId];
	    if(midId-start <= m_min_cell_size || level >= PCA_KDTREE_MAX_DEPTH)
	    {
	        child.leaf = 1;
	        child.start = start;
	        child.size = midId-start;
	    }
	    else
	    {
	        child.leaf = 0;
	        this->build_rec(childId, start, midId, level+1);
	    }
	}
	{
	    // right child
	    int childId = nodes[node_id].firstChildId+1;
	    KdTreeNode& child = nodes[childId];
	    if(end-midId <= m_min_cell_size || level >= PCA_KDTREE_MAX_DEPTH)
	    {
	        child.leaf = 1;
	        child.start = midId;
	        child.size = end-midId;
	    }
	    else
	    {
	        child.leaf = 0;
	        this->build_rec(childId, midId, end, level+1);
	    }
	}
}

template<class DataPoint>
inline int Ponca::KdTree<DataPoint>::partition(int start, int end, int dim, Scalar value)
{
	const auto& points = *m_points.get();
	auto& indices  = *m_indices.get();
	
	auto it = std::partition(indices.begin()+start, indices.begin()+end, [&](int i)
	{
	    return points[i][dim] < value;
	});
	    
	auto distance = std::distance(m_indices->begin(), it);
	
	return static_cast<int>(distance);
}

template<class DataPoint>
inline Ponca::KdTreeRangeIndexQuery<DataPoint::Scalar> Ponca::KdTree<DataPoint>::range_neighbors(int index, DataPoint::Scalar r) const
{
	return RangeIndexQuery(r, index);
}

//
//// Query -----------------------------------------------------------------------
//
//KdTreeKNearestPointQuery KdTree::k_nearest_neighbors(const Vector3& point, int k) const
//{
//	return KNearestPointQuery(this, k, point);
//}
//
//KdTreeKNearestIndexQuery KdTree::k_nearest_neighbors(int index, int k) const
//{
//	return KNearestIndexQuery(this, k, index);
//}
//
//KdTreeNearestPointQuery KdTree::nearest_neighbor(const Vector3& point) const
//{
//	return NearestPointQuery(this, point);
//}
//
//KdTreeNearestIndexQuery KdTree::nearest_neighbor(int index) const
//{
//	return NearestIndexQuery(this, index);
//}
//
//KdTreeRangePointQuery KdTree::range_neighbors(const Vector3& point, Scalar r) const
//{
//	return RangePointQuery(this, r, point);
//}
//
//KdTreeRangeIndexQuery KdTree::range_neighbors(int index, Scalar r) const
//{
//	return RangeIndexQuery(this, r, index);
//}
//
//// Empty Query -----------------------------------------------------------------
//
//KdTreeKNearestPointQuery KdTree::k_nearest_point_query(int k) const
//{
//	return KNearestPointQuery(this, k);
//}
//
//KdTreeKNearestIndexQuery KdTree::k_nearest_index_query(int k) const
//{
//	return KNearestIndexQuery(this, k);
//}
//
//KdTreeNearestPointQuery KdTree::nearest_point_query() const
//{
//	return NearestPointQuery(this);
//}
//
//KdTreeNearestIndexQuery KdTree::nearest_index_query() const
//{
//	return NearestIndexQuery(this);
//}
//
//KdTreeRangePointQuery KdTree::range_point_query(Scalar r) const
//{
//	return RangePointQuery(this, r);
//}
//
//KdTreeRangeIndexQuery KdTree::range_index_query(Scalar r) const
//{
//	return RangeIndexQuery(this, r);
//}
