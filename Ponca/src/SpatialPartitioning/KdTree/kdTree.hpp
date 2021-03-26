/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// KdTree ----------------------------------------------------------------------

template<class DataPoint>
int KdTree<DataPoint>::node_count() const
{
	return static_cast<int>(m_nodes->size());
}

template<class DataPoint>
int KdTree<DataPoint>::index_count() const
{
	return m_indices.size();
}

template<class DataPoint>
int KdTree<DataPoint>::point_count() const
{
	return static_cast<int>(m_points.size());
}

template<class DataPoint>
void KdTree<DataPoint>::clear()
{
	m_points.clear();
	m_nodes.clear();
	m_indices.clear();
}

//template<class DataPoint>
//void KdTree<DataPoint>::build(std::shared_ptr<Vector>& points)
//{
//	std::vector<int> ids;
//	iota(ids.begin(), ids.end(), 0);
//	this->build( points, ids);
//}
//
//template<class DataPoint>
//void KdTree<DataPoint>::build(std::shared_ptr<Vector>& points, const std::vector<int>& sampling)
//{
//	this->clear();
//	
//	m_points = points;
//	
//	m_nodes = std::make_shared<std::vector<KdTreeNode<Scalar>>>();
//	m_nodes->reserve(4 * m_points->size() / m_min_cell_size);
//	m_nodes->emplace_back();
//	m_nodes->back().leaf = false;
//	
//	m_indices = std::make_shared<std::vector<int>>(sampling);//move operator ou std copy
//	
//	int end = static_cast<int>(m_indices->size());
//	
//	this->build_rec(0, 0, end, 1);
//	
////    PCA_DEBUG_ASSERT(this->valid());
//}
//
//template<class DataPoint>
//void KdTree<DataPoint>::rebuild(const std::vector<int>& sampling)
//{
//	//    PCA_DEBUG_ASSERT(sampling.size() <= m_points->size());
//
//    m_nodes->clear();
//    m_nodes->emplace_back();
//    m_nodes->back().leaf = false;
//
//    *m_indices = sampling;
//
//    int end = static_cast<int>(m_indices->size());
//    this->build_rec(0, 0, end, 1);
//
////    PCA_DEBUG_ASSERT(this->valid());
//}


template<class DataPoint>
template<typename VectorUserContainer>
inline void KdTree<DataPoint>::build(const VectorUserContainer & points)
{
	std::vector<int> ids;
	iota(ids.begin(), ids.end(), 0);
	this->build(points, ids);
}

template<class DataPoint>
template<typename VectorUserContainer, typename IndexUserContainer>
inline void KdTree<DataPoint>::build(const VectorUserContainer & points, const IndexUserContainer & sampling)
{
	this->clear();

	//m_points = points;
	std::copy(points.begin(), points.end(), m_points.begin());

	m_nodes = NodeContainer();
	m_nodes.reserve(4 * m_points.size() / m_min_cell_size);
	m_nodes.emplace_back();
	m_nodes.back().leaf = false;

	m_indices = IndexContainer(sampling);//move operator ou std copy

	int end = m_indices.size();

	this->build_rec(0, 0, end, 1);

	   // PCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
template<typename IndexUserContainer>
inline void KdTree<DataPoint>::rebuild(const IndexUserContainer & sampling)
{
	//    PCA_DEBUG_ASSERT(sampling.size() <= m_points->size());

	m_nodes->clear();
	m_nodes->emplace_back();
	m_nodes->back().leaf = false;

	m_indices = sampling;

	int end = static_cast<int>(m_indices.size());
	this->build_rec(0, 0, end, 1);

	//    PCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
bool KdTree<DataPoint>::valid() const
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
		const KdTreeNode<Scalar>& node = m_nodes->operator [](n);
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
std::string KdTree<DataPoint>::to_string() const
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
	    const KdTreeNode<Scalar>& node = m_nodes->operator[](n);
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
int KdTree<DataPoint>::min_cell_size() const
{
	return m_min_cell_size;
}

template<class DataPoint>
void KdTree<DataPoint>::set_min_cell_size(int min_cell_size)
{
	m_min_cell_size = min_cell_size;
}

template<class DataPoint>
void KdTree<DataPoint>::build_rec(int node_id, int start, int end, int level)
{
	auto& nodes = m_nodes;
	const auto& points  = m_points;
	const auto& indices = m_indices;
	
	KdTreeNode<Scalar>& node = nodes[node_id];
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
	    KdTreeNode<Scalar> n;
	    n.size = 0;
	    nodes.push_back(n);
	    nodes.push_back(n);
	}
	{
	    // left child
	    int childId = nodes[node_id].firstChildId;
	    KdTreeNode<Scalar>& child = nodes[childId];
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
	    KdTreeNode<Scalar>& child = nodes[childId];
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
int KdTree<DataPoint>::partition(int start, int end, int dim, Scalar value)
{
	const auto& points = m_points;
	auto& indices  = m_indices;
	
	auto it = std::partition(indices.begin()+start, indices.begin()+end, [&](int i)
	{
	    return points[i][dim] < value;
	});
	    
	auto distance = std::distance(m_indices.begin(), it);
	
	return static_cast<int>(distance);
}

