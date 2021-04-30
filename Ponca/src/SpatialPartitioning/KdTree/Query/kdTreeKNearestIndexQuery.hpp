/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint>
KdTreeKNearestIndexIterator<DataPoint> KdTreeKNearestIndexQuery<DataPoint>::begin()
{
    this->search();
    return KdTreeKNearestIndexIterator<DataPoint>(QueryType::m_queue.begin());
}

template<class DataPoint>
KdTreeKNearestIndexIterator<DataPoint> KdTreeKNearestIndexQuery<DataPoint>::end()
{
    return KdTreeKNearestIndexIterator<DataPoint>(QueryType::m_queue.end());
}

template<class DataPoint>
void KdTreeKNearestIndexQuery<DataPoint>::search()
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = points[QueryType::m_index];

    QueryAccelType::m_stack.clear();
    QueryAccelType::m_stack.push({0,0});

    QueryType::m_queue.clear();
    QueryType::m_queue.push({-1,std::numeric_limits<Scalar>::max()});

    while(!QueryAccelType::m_stack.empty())
    {
        auto& qnode = QueryAccelType::m_stack.top();
        const auto& node  = nodes[qnode.index];

        if(qnode.squared_distance < QueryType::m_queue.bottom().squared_distance)
        {
            if(node.leaf)
            {
                QueryAccelType::m_stack.pop();
                int end = node.start + node.size;
                for(int i=node.start; i<end; ++i)
                {
                    int idx = indices[i];
                    if(QueryType::m_index == idx) continue;

                    Scalar d = (point - points[idx]).squaredNorm();
                    QueryType::m_queue.push({idx, d});
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
                qnode.squared_distance         = newOff * newOff;
            }
        }
        else
        {
            QueryAccelType::m_stack.pop();
        }
    }
}