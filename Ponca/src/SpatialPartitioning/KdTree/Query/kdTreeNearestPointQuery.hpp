/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <typename Traits>
auto KdTreeNearestPointQuery<Traits>::begin()
    -> KdTreeNearestIterator<IndexType>
{
    QueryAccelType::reset();
    QueryType::reset();
    this->search();
    return KdTreeNearestIterator<IndexType>(QueryType::m_nearest);
}

template <typename Traits>
auto KdTreeNearestPointQuery<Traits>::end()
    -> KdTreeNearestIterator<IndexType>
{
    return KdTreeNearestIterator<IndexType>(QueryType::m_nearest + 1);
}

template <typename Traits>
void KdTreeNearestPointQuery<Traits>::search()
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = QueryType::input();

    if (nodes.empty() || points.empty() || indices.empty())
        throw std::invalid_argument("Empty KdTree");

    while(!QueryAccelType::m_stack.empty())
    {
        auto& qnode = QueryAccelType::m_stack.top();
        const auto& node  = nodes[qnode.index];

        if(qnode.squared_distance < QueryType::m_squared_distance)
        {
            if(node.is_leaf())
            {
                QueryAccelType::m_stack.pop();
                IndexType end = node.leaf.start + node.leaf.size;
                for(IndexType i=node.leaf.start; i<end; ++i)
                {
                    IndexType idx = indices[i];
                    Scalar d = (point - points[idx].pos()).squaredNorm();
                    if(d < QueryType::m_squared_distance)
                    {
                        QueryType::m_nearest = idx;
                        QueryType::m_squared_distance = d;
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
}
