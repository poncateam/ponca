/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"
#include "./KdTreeNearestPointQuery.h"

class KdTree;

namespace Ponca {

template <typename VectorType>
KdTreeNearestPointQuery<VectorType>::KdTreeNearestPointQuery() :
    KdTreeQuery(),
    NearestPointQuery<VectorType>()
{
    cout << "Test" << endl;
}

template <typename VectorType>
KdTreeNearestPointQuery<VectorType>::KdTreeNearestPointQuery(const KdTree* kdtree) :
    KdTreeQuery(kdtree),
    NearestPointQuery<VectorType>()
{
}

template <typename VectorType>
KdTreeNearestPointQuery<VectorType>::KdTreeNearestPointQuery(const KdTree* kdtree, const VectorType& point) :
    KdTreeQuery(kdtree),
    NearestPointQuery<VectorType>(point)
{
}

template <typename VectorType>
KdTreeNearestPointIterator KdTreeNearestPointQuery<VectorType>::begin()
{
    this->search();
    return KdTreeNearestPointIterator(m_nearest);
}

template <typename VectorType>
KdTreeNearestPointIterator KdTreeNearestPointQuery<VectorType>::end()
{
    return KdTreeNearestPointIterator(m_nearest+1);
}

template <typename VectorType>
void KdTreeNearestPointQuery<VectorType>::search()
{
    const auto& nodes   = m_kdtree->node_data();
    const auto& points  = m_kdtree->point_data();
    const auto& indices = m_kdtree->index_data();

    if (nodes.empty() || points.empty() || indices.empty())
        throw invalid_argument("Empty KdTree");

    m_stack.clear();
    m_stack.push({0,0});

    m_nearest = indices[0];
    m_squared_distance = (m_point - points[m_nearest]).squaredNorm();

    while(!m_stack.empty())
    {
        auto& qnode = m_stack.top();
        const auto& node  = nodes[qnode.index];

        if(qnode.squared_distance < m_squared_distance)
        {
            if(node.leaf)
            {
                m_stack.pop();
                int end = node.start + node.size;
                for(int i=node.start; i<end; ++i)
                {
                    int idx = indices[i];

                    Scalar d = (m_point - points[idx]).squaredNorm();
                    if(d < m_squared_distance)
                    {
                        m_nearest = idx;
                        m_squared_distance = d;
                    }
                }
            }
            else
            {
                // replace the stack top by the farthest and push the closest
                Scalar newOff = m_point[node.dim] - node.splitValue;
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
}

} // namespace ponca
