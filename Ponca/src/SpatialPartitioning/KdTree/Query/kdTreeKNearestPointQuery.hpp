/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"
#include "./KdTreeKNearestPointQuery.h"

class KdTree;

namespace Ponca {

template <typename VectorType>
KdTreeKNearestPointQuery<VectorType>::KdTreeKNearestPointQuery() :
    KdTreeQuery(),
    KNearestPointQuery<VectorType>()
{
    cout << "Test" << endl;
}

template <typename VectorType>
KdTreeKNearestPointQuery<VectorType>::KdTreeKNearestPointQuery(const KdTree* kdtree, int k) :
    KdTreeQuery(kdtree),
    KNearestPointQuery<VectorType>(k)
{
}

template <typename VectorType>
KdTreeKNearestPointQuery<VectorType>::KdTreeKNearestPointQuery(const KdTree* kdtree, int k, const VectorType& point) :
    KdTreeQuery(kdtree),
    KNearestPointQuery<VectorType>(k, point)
{
}

template <typename VectorType>
KdTreeKNearestPointIterator KdTreeKNearestPointQuery<VectorType>::begin()
{
    this->search();
    return KdTreeKNearestPointIterator(KNearestQuery::m_queue.begin());
}

template <typename VectorType>
KdTreeKNearestPointIterator KdTreeKNearestPointQuery<VectorType>::end()
{
    return KdTreeKNearestPointIterator(KNearestQuery::m_queue.end());
}

template <typename VectorType>
void KdTreeKNearestPointQuery<VectorType>::search()
{
    const auto& nodes   = m_kdtree->node_data();
    const auto& points  = m_kdtree->point_data();
    const auto& indices = m_kdtree->index_data();

    m_stack.clear();
    m_stack.push({0,0});

    KNearestQuery::m_queue.clear();
    KNearestQuery::m_queue.push({-1,std::numeric_limits<Scalar>::max()});

    while(!m_stack.empty())
    {
        auto& qnode = m_stack.top();
        const auto& node  = nodes[qnode.index];

        if(qnode.squared_distance < KNearestQuery::m_queue.bottom().squared_distance)
        {
            if(node.leaf)
            {
                m_stack.pop();
                int end = node.start + node.size;
                for(int i=node.start; i<end; ++i)
                {
                    int idx = indices[i];

                    Scalar d = (PointQuery<VectorType>::m_point - points[idx]).squaredNorm();
                    KNearestQuery::m_queue.push({idx, d});
                }
            }
            else
            {
                // replace the stack top by the farthest and push the closest
                Scalar newOff = PointQuery<VectorType>::m_point[node.dim] - node.splitValue;
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

}