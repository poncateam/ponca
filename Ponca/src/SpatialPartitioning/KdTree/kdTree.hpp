/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// KdTree ----------------------------------------------------------------------

template<typename Traits>
template<typename PointUserContainer, typename PointConverter>
void KdTreeBase<Traits>::build(PointUserContainer&& points, PointConverter c)
{
    IndexContainer ids (points.size());
    std::iota(ids.begin(), ids.end(), IndexType(0));
    this->buildWithSampling(std::forward<PointUserContainer>(points), std::move(ids), std::move(c));
}

template<typename Traits>
bool StaticKdTreeBase<Traits>::valid() const
{
    if (m_bufs.points_size == 0)
        return m_bufs.nodes_size == 0 && m_bufs.indices_size == 0;

    if(m_bufs.nodes_size == 0 || m_bufs.indices_size == 0)
    {
        std::cerr << "KdTree validation check failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
        return false;
    }

    std::vector<bool> b(pointCount(), false);
    for(unsigned int i = 0; i < sampleCount(); ++i)
    {
        const int idx = m_bufs.indices[i];
        if(idx < 0 || pointCount() <= idx || b[idx])
        {
            std::cerr << "KdTree validation check failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
            return false;
        }
        b[idx] = true;
    }

    for(NodeIndexType n=0;n<nodeCount();++n)
    {
        const NodeType& node = m_bufs.nodes[n];
        if(node.is_leaf())
        {
            if(sampleCount() <= node.leaf_start() || node.leaf_start()+node.leaf_size() > sampleCount())
            {
                std::cerr << "KdTree validation check failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
                return false;
            }
        }
        else
        {
            if(node.inner_split_dim() < 0 || DataPoint::Dim-1 < node.inner_split_dim())
            {
                std::cerr << "KdTree validation check failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
                return false;
            }
            if(nodeCount() <= node.inner_first_child_id() || nodeCount() <= node.inner_first_child_id()+1)
            {
                std::cerr << "KdTree validation check failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
                return false;
            }
        }
    }

    return true;
}

template<typename Traits>
void StaticKdTreeBase<Traits>::print(std::ostream& os, bool verbose) const
{
    os << "KdTree:";
    os << "\n  MaxNodes: " << MAX_NODE_COUNT;
    os << "\n  MaxPoints: " << MAX_POINT_COUNT;
    os << "\n  MaxDepth: " << MAX_DEPTH;
    os << "\n  PointCount: " << pointCount();
    os << "\n  SampleCount: " << sampleCount();
    os << "\n  NodeCount: " << nodeCount();

    if (!verbose)
    {
        return;
    }

    os << "\n  Samples: [";
    static constexpr IndexType SAMPLES_PER_LINE = 10;
    for (IndexType i = 0; i < sampleCount(); ++i)
    {
        os << (i == 0 ? "" : ",");
        os << (i % SAMPLES_PER_LINE == 0 ? "\n    " : " ");
        os << m_bufs.indices[i];
    }

    os << "]\n  Nodes:";
    for (NodeIndexType n = 0; n < nodeCount(); ++n)
    {
        const NodeType& node = m_bufs.nodes[n];
        if (node.is_leaf())
        {
            os << "\n    - Type: Leaf";
            os << "\n      Start: " << node.leaf_start();
            os << "\n      Size: " << node.leaf_size();
        }
        else
        {
            os << "\n    - Type: Inner";
            os << "\n      SplitDim: " << node.inner_split_dim();
            os << "\n      SplitValue: " << node.inner_split_value();
            os << "\n      FirstChild: " << node.inner_first_child_id();
        }
    }
}

template<typename Traits>
template<typename PointUserContainer, typename IndexUserContainer, typename PointConverter>
void KdTreeBase<Traits>::buildWithSampling(
    PointUserContainer&& points, IndexUserContainer&& sampling, PointConverter c
) {
    PONCA_DEBUG_ASSERT(static_cast<IndexType>(pointCount()) <= MAX_POINT_COUNT);
    Base::m_leaf_count = 0;

    // Move, copy or convert input samples
    c(std::forward<PointUserContainer>(points), Base::m_bufs.points);
    Base::m_bufs.points_size = points.size();

    Base::m_bufs.indices_size = sampling.size();
    Base::m_bufs.indices      = std::move(sampling);

    Base::m_bufs.nodes.reserve(4 * Base::pointCount() / Base::m_min_cell_size);
    Base::m_bufs.nodes.emplace_back();

    this->buildRec(0, 0, Base::sampleCount(), 1);
    Base::m_bufs.nodes_size   = Base::m_bufs.nodes.size();

    PONCA_DEBUG_ASSERT(this->valid());
}

template<typename Traits>
void KdTreeBase<Traits>::buildRec(NodeIndexType node_id, IndexType start, IndexType end, int level)
{
    NodeType& node = Base::m_bufs.nodes[node_id];
    AabbType aabb;
    for(IndexType i=start; i<end; ++i)
        aabb.extend(Base::m_bufs.points[Base::m_bufs.indices[i]].pos());

    node.set_is_leaf(
        end-start <= Base::m_min_cell_size ||
        level >= Traits::MAX_DEPTH ||
        // Since we add 2 nodes per inner node we need to stop if we can't add
        // them both
        static_cast<NodeIndexType>(Base::m_bufs.nodes.size()) > Base::MAX_NODE_COUNT - 2);

    node.configure_range(start, end-start, aabb);
    if (node.is_leaf())
    {
        ++Base::m_leaf_count;
    }
    else
    {
        IndexType split_dim = 0;
        (Scalar(0.5) * aabb.diagonal()).maxCoeff(&split_dim);
        node.configure_inner(aabb.center()[split_dim], static_cast<IndexType>(Base::m_bufs.nodes.size()), split_dim);
        Base::m_bufs.nodes.emplace_back();
        Base::m_bufs.nodes.emplace_back();

        IndexType mid_id = this->partition(start, end, split_dim, node.inner_split_value());
        buildRec(node.inner_first_child_id(), start, mid_id, level+1);
        buildRec(node.inner_first_child_id()+1, mid_id, end, level+1);
    }
}

template<typename Traits>
auto KdTreeBase<Traits>::partition(IndexType start, IndexType end, int dim, Scalar value)
    -> IndexType
{
    const auto& points = Base::m_bufs.points;
    
    auto it = std::partition(std::begin(Base::m_bufs.indices)+start, std::begin(Base::m_bufs.indices)+end, [&](IndexType i)
    {
        return points[i].pos()[dim] < value;
    });

    auto distance = std::distance(std::begin(Base::m_bufs.indices), it);
    
    return static_cast<IndexType>(distance);
}
