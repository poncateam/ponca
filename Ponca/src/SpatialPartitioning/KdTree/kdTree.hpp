/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// KdTree ----------------------------------------------------------------------

template<class DataPoint>
int KdTree<DataPoint>::node_count() const
{
    return static_cast<int>(m_nodes.size());
}

template<class DataPoint>
int KdTree<DataPoint>::index_count() const
{
    return static_cast<int>(m_indices.size());
}

template<class DataPoint>
int KdTree<DataPoint>::point_count() const
{
    return static_cast<int>(m_points.size());
}


template<class DataPoint>
int KdTree<DataPoint>::leaf_count() const
{
    return m_leaf_count;
}

template<class DataPoint>
void KdTree<DataPoint>::clear()
{
    m_points.clear();
    m_nodes.clear();
    m_indices.clear();
    m_leaf_count = 0;
}

template<class DataPoint>
template<typename PointUserContainer, typename Converter>
inline void KdTree<DataPoint>::build(const PointUserContainer& points, Converter c)
{
    std::vector<int> ids(points.size());
    std::iota(ids.begin(), ids.end(), 0);
    this->buildWithSampling(points, ids, c);
}

template<class DataPoint>
template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
inline void KdTree<DataPoint>::buildWithSampling(const PointUserContainer& points,
                                                 const IndexUserContainer& sampling,
                                                 Converter c)
{
    this->clear();

    // Copy or convert input samples
    c( points, m_points );

    m_nodes = NodeContainer();
    m_nodes.reserve(4 * point_count() / m_min_cell_size);
    m_nodes.emplace_back();
    m_nodes.back().leaf = false;

    m_indices = IndexContainer(sampling);//move operator ou std copy

    this->build_rec(0, 0, index_count(), 1);

    PONCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
template<typename IndexUserContainer>
inline void KdTree<DataPoint>::rebuild(const IndexUserContainer & sampling)
{
    PONCA_DEBUG_ASSERT(sampling.size() <= m_points->size());

    m_nodes.clear();
    m_nodes.emplace_back();
    m_nodes.back().leaf = false;

    m_indices = sampling;

    this->build_rec(0, 0, index_count(), 1);

    PONCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint>
bool KdTree<DataPoint>::valid() const
{
    PONCA_DEBUG_ERROR;
    return false;

    if (m_points.empty())
        return m_nodes.empty() && m_indices.empty();
        
    if(m_nodes.empty() || m_indices.empty())
    {
        PONCA_DEBUG_ERROR;
        return false;
    }
        
    if(point_count() < index_count())
    {
        PONCA_DEBUG_ERROR;
        return false;
    }
        
    std::vector<bool> b(point_count(), false);
    for(int idx : m_indices)
    {
        if(idx < 0 || point_count() <= idx || b[idx])
        {
            PONCA_DEBUG_ERROR;
            return false;
        }
        b[idx] = true;
    }
        
    for(int n=0;n<node_count();++n)
    {
        const KdTreeNode<Scalar>& node = m_nodes.operator[](n);
        if(node.leaf)
        {
            if(index_count() <= node.start || index_count() < node.start+node.size)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
        }
        else
        {
            if(node.dim < 0 || 2 < node.dim)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
            if(node_count() <= node.firstChildId || node_count() <= node.firstChildId+1u)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
        }
    }

    return true;
}

template<class DataPoint>
std::string KdTree<DataPoint>::to_string() const
{
    if (m_indices.empty()) return "";
    
    std::stringstream str;
    str << "indices (" << index_count() << ") :\n";
    for(int i=0; i<index_count(); ++i)
    {
        str << "  " << i << ": " << m_indices.operator[](i) << "\n";
    }
    str << "nodes (" << node_count() << ") :\n";
    for(int n=0; n< node_count(); ++n)
    {
        const KdTreeNode<Scalar>& node = m_nodes.operator[](n);
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
    KdTreeNode<Scalar>& node = m_nodes[node_id];
    Aabb aabb;
    for(int i=start; i<end; ++i)
        aabb.extend(m_points[m_indices[i]].pos());
    
    int dim;
    (Scalar(0.5) * (aabb.max() - aabb.min())).maxCoeff(&dim);
    node.dim = dim;
    node.splitValue = aabb.center()(dim);
    
    int midId = this->partition(start, end, dim, node.splitValue);
    node.firstChildId = m_nodes.size();
    
    {
        KdTreeNode<Scalar> n;
        n.size = 0;
        m_nodes.push_back(n);
        m_nodes.push_back(n);
    }
    {
        // left child
        int childId = m_nodes[node_id].firstChildId;
        KdTreeNode<Scalar>& child = m_nodes[childId];
        if(midId-start <= m_min_cell_size || level >= PCA_KDTREE_MAX_DEPTH)
        {
            child.leaf = 1;
            child.start = start;
            child.size = midId-start;
            m_leaf_count++;
        }
        else
        {
            child.leaf = 0;
            this->build_rec(childId, start, midId, level+1);
        }
    }
    {
        // right child
        int childId = m_nodes[node_id].firstChildId+1;
        KdTreeNode<Scalar>& child = m_nodes[childId];
        if(end-midId <= m_min_cell_size || level >= PCA_KDTREE_MAX_DEPTH)
        {
            child.leaf = 1;
            child.start = midId;
            child.size = end-midId;
            m_leaf_count++;
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
        return points[i].pos()[dim] < value;
    });
        
    auto distance = std::distance(m_indices.begin(), it);
    
    return static_cast<int>(distance);
}

