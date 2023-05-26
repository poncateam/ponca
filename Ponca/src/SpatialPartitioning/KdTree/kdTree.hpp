/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// KdTree ----------------------------------------------------------------------

template<class DataPoint, class Compatibility>
void KdTree<DataPoint, Compatibility>::clear()
{
    m_points.clear();
    m_nodes.clear();
    m_indices.clear();
    m_leaf_count = 0;
}

template<class DataPoint, class Compatibility>
template<typename PointUserContainer, typename Converter>
inline void KdTree<DataPoint, Compatibility>::build(const PointUserContainer& points, Converter c)
{
    IndexContainer ids(points.size());
    std::iota(ids.begin(), ids.end(), 0);
    this->buildWithSampling(points, std::move(ids), c);
}

template<class DataPoint, class Compatibility>
template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
inline void KdTree<DataPoint, Compatibility>::buildWithSampling(const PointUserContainer& points,
                                                                const IndexUserContainer& sampling,
                                                                Converter c)
{
    this->clear();

    // Copy or convert input samples
    c( points, m_points );

    m_nodes = NodeContainer();
    m_nodes.reserve(4 * point_count() / m_min_cell_size);
    m_nodes.emplace_back();

    m_indices = IndexContainer(sampling);//move operator ou std copy

    this->build_rec(0, 0, index_count(), 1);

    PONCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint, class Compatibility>
template<typename IndexUserContainer>
inline void KdTree<DataPoint, Compatibility>::rebuild(const IndexUserContainer & sampling)
{
    PONCA_DEBUG_ASSERT(sampling.size() <= m_points->size());

    m_nodes.clear();
    m_nodes.emplace_back();

    m_indices = sampling;

    this->build_rec(0, 0, index_count(), 1);

    PONCA_DEBUG_ASSERT(this->valid());
}

template<class DataPoint, class Compatibility>
bool KdTree<DataPoint, Compatibility>::valid() const
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

    for(NodeCountType n=0;n<node_count();++n)
    {
        const NodeType& node = m_nodes.operator[](n);
        if(node.is_leaf())
        {
            if(index_count() <= node.leaf.start || index_count() < node.leaf.start+node.leaf.size)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
        }
        else
        {
            if(node.inner.dim < 0 || 2 < node.inner.dim)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
            if(node_count() <= node.inner.first_child_id || node_count() <= node.inner.first_child_id+1)
            {
                PONCA_DEBUG_ERROR;
                return false;
            }
        }
    }

    return true;
}

template<class DataPoint, class Compatibility>
std::string KdTree<DataPoint, Compatibility>::to_string() const
{
    if (m_indices.empty()) return "";
    
    std::stringstream str;
    str << "indices (" << index_count() << ") :\n";
    for(IndexCountType i=0; i<index_count(); ++i)
    {
        str << "  " << i << ": " << m_indices.operator[](i) << "\n";
    }
    str << "nodes (" << node_count() << ") :\n";
    for(NodeCountType n=0; n< node_count(); ++n)
    {
        const NodeType& node = m_nodes.operator[](n);
        if(node.is_leaf())
        {
            auto end = node.leaf.start + node.leaf.size;
            str << "  leaf: start=" << node.leaf.start << " end=" << end << " (size=" << node.leaf.size << ")\n";
        }
        else
        {
            str << "  node: dim=" << node.inner.dim << " split=" << node.inner.split_value << " child=" << node.inner.first_child_id << "\n";
        }
    }
    return str.str();
}

template<class DataPoint, class Compatibility>
void KdTree<DataPoint, Compatibility>::build_rec(NodeCountType node_id, IndexCountType start, IndexCountType end, DepthType level)
{
    NodeType& node = m_nodes[node_id];
    for(IndexCountType i=start; i<end; ++i)
        node.aabb.extend(m_points[m_indices[i]].pos());

    if (end-start <= m_min_cell_size || level >= PCA_KDTREE_MAX_DEPTH)
    {
        PONCA_ASSERT_MSG(end-start <= std::numeric_limits<LeafSizeType>::max(), "Leaf size overflow");

        node.set_is_leaf(true);

        node.leaf.start = start;
        node.leaf.size = static_cast<LeafSizeType>(end-start);
        ++m_leaf_count;
    }
    else
    {
        node.set_is_leaf(false);

        DimType dim;
        (Scalar(0.5) * (node.aabb.max() - node.aabb.min())).maxCoeff(&dim);
        node.inner.dim = dim;
        node.inner.split_value = node.aabb.center()(dim);

        IndexCountType mid_id = this->partition(start, end, dim, node.inner.split_value);
        node.inner.first_child_id = m_nodes.size();
        m_nodes.emplace_back();
        m_nodes.emplace_back();

        build_rec(node.inner.first_child_id, start, mid_id, level+1);
        build_rec(node.inner.first_child_id+1, mid_id, end, level+1);
    }
}

template<class DataPoint, class Compatibility>
auto KdTree<DataPoint, Compatibility>::partition(IndexCountType start, IndexCountType end, DimType dim, Scalar value)
    -> IndexCountType
{
    const auto& points = m_points;
    auto& indices  = m_indices;
    
    auto it = std::partition(indices.begin()+start, indices.begin()+end, [&](IndexType i)
    {
        return points[i].pos()[dim] < value;
    });
        
    auto distance = std::distance(m_indices.begin(), it);
    
    return static_cast<IndexCountType>(distance);
}
