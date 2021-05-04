/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint>
KdTreeNearestPointIterator KdTreeNearestPointQuery<DataPoint>::begin()
{
    this->search();
    return KdTreeNearestPointIterator( QueryType::m_nearest );
}

template <class DataPoint>
KdTreeNearestPointIterator KdTreeNearestPointQuery<DataPoint>::end()
{
    return KdTreeNearestPointIterator( QueryType::m_nearest + 1);
}

template <class DataPoint>
void KdTreeNearestPointQuery<DataPoint>::search()
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();

    if (nodes.empty() || points.empty() || indices.empty())
        throw std::invalid_argument("Empty KdTree");

    QueryAccelType::m_stack.clear();
    QueryAccelType::m_stack.push({0,0});

    QueryType::m_nearest = indices[0];
    QueryType::m_squared_distance = (QueryType::m_point.pos() - points[QueryType::m_nearest].pos()).squaredNorm();

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
                    Scalar d = (QueryType::m_point.pos() - points[idx].pos()).squaredNorm();
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
                Scalar newOff = QueryType::m_point.pos()[node.dim] - node.splitValue;
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
