/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint, class Adapter>
KdTreeNearestIterator KdTreeNearestIndexQuery<DataPoint, Adapter>::begin()
{
    QueryAccelType::reset();
    QueryType::reset();
    this->search();
    return KdTreeNearestIterator(QueryType::m_nearest);
}

template <class DataPoint, class Adapter>
KdTreeNearestIterator KdTreeNearestIndexQuery<DataPoint, Adapter>::end()
{
    return KdTreeNearestIterator(QueryType::m_nearest + 1);
}

template <class DataPoint, class Adapter>
void KdTreeNearestIndexQuery<DataPoint, Adapter>::search()
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = points[QueryType::input()].pos();

    if (nodes.empty() || points.empty() || indices.empty())
        throw std::invalid_argument("Empty KdTree");

    while(!QueryAccelType::m_stack.empty())
    {
        auto& qnode = QueryAccelType::m_stack.top();
        const auto& node  = nodes[qnode.index];

        if(qnode.squared_distance < QueryType::m_squared_distance)
        {
            if(node.leaf)
            {
                QueryAccelType::m_stack.pop();
                int end = node.start + node.size;
                for(int i=node.start; i<end; ++i)
                {
                    int idx = indices[i];
                    if(QueryType::input() == idx) continue;

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
                Scalar newOff = point[node.dim] - node.splitValue;
                QueryAccelType::m_stack.push();
                if(newOff < 0)
                {
                    QueryAccelType::m_stack.top().index = node.firstChildId;
                    qnode.index         = node.firstChildId+1;
                }
                else
                {
                    QueryAccelType::m_stack.top().index = node.firstChildId+1;
                    qnode.index         = node.firstChildId;
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
