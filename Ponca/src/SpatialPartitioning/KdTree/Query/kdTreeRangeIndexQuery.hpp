/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint>
KdTreeRangeIndexIterator<DataPoint> KdTreeRangeIndexQuery<DataPoint>::begin()
{
    KdTreeRangeIndexIterator<DataPoint> it(this);
    this->initialize(it);
    this->advance(it);
    return it;
}

template<class DataPoint>
KdTreeRangeIndexIterator<DataPoint> KdTreeRangeIndexQuery<DataPoint>::end()
{
    return KdTreeRangeIndexIterator<DataPoint>(this, m_kdtree->point_count());
}

template<class DataPoint>
void KdTreeRangeIndexQuery<DataPoint>::initialize(KdTreeRangeIndexIterator<DataPoint>& it)
{
    m_stack.clear();
    m_stack.push();
    m_stack.top().index = 0;
    m_stack.top().squared_distance = 0;
    it.m_index = -1;
    it.m_start = 0;
    it.m_end   = 0;
}

template<class DataPoint>
void KdTreeRangeIndexQuery<DataPoint>::advance(KdTreeRangeIndexIterator<DataPoint>& it)
{
    const auto& nodes   = m_kdtree->node_data();
    const auto& points  = m_kdtree->point_data();
    const auto& indices = m_kdtree->index_data();
    const auto& point   = points[m_index];

    for(int i=it.m_start; i<it.m_end; ++i)
    {
        int idx = indices[i];
        if(idx == m_index) continue;

        Scalar d = (point - points[idx]).squaredNorm();
        if(d < m_squared_radius)
        {
            it.m_index = idx;
            it.m_start = i+1;
            return;
        }
    }

    while(!m_stack.empty())
    {
        auto& qnode = m_stack.top();
        const auto& node = nodes[qnode.index];

        if(qnode.squared_distance < m_squared_radius)
        {
            if(node.leaf)
            {
                m_stack.pop();
                it.m_start = node.start;
                it.m_end   = node.start + node.size;
                for(int i=it.m_start; i<it.m_end; ++i)
                {
                    int idx = indices[i];
                    if(idx == m_index) continue;

                    Scalar d = (point - points[idx]).squaredNorm();
                    if(d < m_squared_radius)
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
                Scalar newOff = point[node.dim] - node.splitValue;
                m_stack.push();
                if(newOff < 0)
                {
                    m_stack.top().index = node.firstChildId;
                    qnode.index         = node.firstChildId+1;
                }
                else
                {
                    m_stack.top().index = node.firstChildId+1;
                    qnode.index         = node.firstChildId;
                }
                m_stack.top().squared_distance = qnode.squared_distance;
                qnode.squared_distance         = newOff*newOff;
            }
        }
        else
        {
            m_stack.pop();
        }
    }
    it.m_index = static_cast<int>(points.size());
}   