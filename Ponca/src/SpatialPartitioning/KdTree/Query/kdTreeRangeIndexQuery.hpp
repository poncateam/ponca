/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint>
typename KdTreeRangeIndexQuery<DataPoint>::Iterator
KdTreeRangeIndexQuery<DataPoint>::begin()
{
    Iterator it(this);
    this->initialize(it);
    this->advance(it);
    return it;
}

template<class DataPoint>
typename KdTreeRangeIndexQuery<DataPoint>::Iterator
KdTreeRangeIndexQuery<DataPoint>::end()
{
    return Iterator(this, QueryAccelType::m_kdtree->point_count());
}

template<class DataPoint>
void KdTreeRangeIndexQuery<DataPoint>::initialize(Iterator& it)
{
    QueryAccelType::m_stack.clear();
    QueryAccelType::m_stack.push();
    QueryAccelType::m_stack.top().index = 0;
    QueryAccelType::m_stack.top().squared_distance = 0;
    it = Iterator(this);
}

template<class DataPoint>
void KdTreeRangeIndexQuery<DataPoint>::advance(Iterator& it)
{
    const auto& nodes   = QueryAccelType::m_kdtree->node_data();
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = points[QueryType::input()];

    for(int i=it.m_start; i<it.m_end; ++i)
    {
        int idx = indices[i];
        if(idx == QueryType::input()) continue;

        Scalar d = (point.pos() - points[idx].pos()).squaredNorm();
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
            if(node.leaf)
            {
                QueryAccelType::m_stack.pop();
                it.m_start = node.start;
                it.m_end   = node.start + node.size;
                for(int i=it.m_start; i<it.m_end; ++i)
                {
                    int idx = indices[i];
                    if(idx == QueryType::input()) continue;

                    Scalar d = (point.pos() - points[idx].pos()).squaredNorm();
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
                Scalar newOff = point.pos()[node.dim] - node.splitValue;
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
    it.m_index = static_cast<int>(points.size());
}   
