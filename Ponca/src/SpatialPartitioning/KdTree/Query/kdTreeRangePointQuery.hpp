/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <typename Traits>
auto KdTreeRangePointQuery<Traits>::begin() -> Iterator
{
    QueryAccelType::reset();
    QueryType::reset();
    Iterator it(this);
    this->advance(it);
    return it;
}

template <typename Traits>
auto KdTreeRangePointQuery<Traits>::end() -> Iterator
{
    return Iterator(this, QueryAccelType::m_kdtree->point_count());
}

template <typename Traits>
void KdTreeRangePointQuery<Traits>::advance(Iterator& it)
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = QueryType::input();

    for(IndexType i=it.m_start; i<it.m_end; ++i)
    {
        IndexType idx = indices[i];

        Scalar d = (point - points[idx].pos()).squaredNorm();
        if(d < QueryType::m_squared_radius)
        {
            it.m_index = idx;
            it.m_start = i+1;
            return;
        }
    }

    while(!QueryAccelType::m_stack.empty())
    {
        auto& qnode = QueryAccelType::m_stack.top();
        const auto& node = nodes[qnode.index];

        if(qnode.squared_distance < QueryType::m_squared_radius)
        {
            if(node.is_leaf())
            {
                QueryAccelType::m_stack.pop();
                it.m_start = node.leaf.start;
                it.m_end   = node.leaf.start + node.leaf.size;
                for(IndexType i=it.m_start; i<it.m_end; ++i)
                {
                    IndexType idx = indices[i];

                    Scalar d = (point - points[idx].pos()).squaredNorm();
                    if(d < QueryType::m_squared_radius)
                    {
                        it.m_index = idx;
                        it.m_start = i+1;
                        return;
                    }
                }
            }
            else
            {
                // replace the stack top by the farthest and push the closest
                Scalar newOff = point[node.inner.dim] - node.inner.split_value;
                QueryAccelType::m_stack.push();
                if(newOff < 0)
                {
                    QueryAccelType::m_stack.top().index = node.inner.first_child_id;
                    qnode.index         = node.inner.first_child_id+1;
                }
                else
                {
                    QueryAccelType::m_stack.top().index = node.inner.first_child_id+1;
                    qnode.index         = node.inner.first_child_id;
                }
                QueryAccelType::m_stack.top().squared_distance = qnode.squared_distance;
                qnode.squared_distance         = newOff*newOff;
            }
        }
        else
        {
            QueryAccelType::m_stack.pop();
        }
    }
    it.m_index = static_cast<IndexType>(points.size());
}
