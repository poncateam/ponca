/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeRangeIterator.h"

namespace Ponca {

template <typename Traits>
class KdTreeRangeIndexQuery : public KdTreeQuery<Traits>,
    public RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = RangeIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeRangeIterator<IndexType, DataPoint, KdTreeRangeIndexQuery>;

protected:
    friend Iterator;

public:

    KdTreeRangeIndexQuery(const KdTreeBase<Traits>* kdtree, Scalar radius, IndexType index) :
        KdTreeQuery<Traits>(kdtree), RangeIndexQuery<IndexType, Scalar>(radius, index)
    {
    }

public:
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        Iterator it(this);
        this->advance(it);
        return it;
    }
    inline Iterator end(){
        return Iterator(this, QueryAccelType::m_kdtree->point_count());
    }

protected:
    inline void advance(Iterator& it){
        const auto& points  = QueryAccelType::m_kdtree->point_data();
        const auto& indices = QueryAccelType::m_kdtree->index_data();
        const auto& point   = QueryAccelType::m_kdtree->point_data()[QueryType::input()].pos();

        auto descentDistanceThreshold = [this](){return QueryType::descentDistanceThreshold();};
        auto skipFunctor              = [this](IndexType idx){return QueryType::input() == idx;};
        auto processNeighborFunctor   = [&it](IndexType idx, IndexType i, Scalar)
        {
            it.m_index = idx;
            it.m_start = i+1;
            return true;
        };

        if (points.empty() || indices.empty())
            throw std::invalid_argument("Empty KdTree");

        for(IndexType i=it.m_start; i<it.m_end; ++i)
        {
            IndexType idx = indices[i];
            if(skipFunctor(idx)) continue;

            Scalar d = (point - points[idx].pos()).squaredNorm();
            if(d < descentDistanceThreshold())
            {
                if( processNeighborFunctor(idx, i, d) ) return;
            }
        }

        if (KdTreeQuery<Traits>::search_internal(point,
                                                 [&it](IndexType start, IndexType end)
                                                 {
                                                     it.m_start = start;
                                                     it.m_end   = end;
                                                 },
                                                 descentDistanceThreshold,
                                                 skipFunctor,
                                                 processNeighborFunctor))
            it.m_index = static_cast<IndexType>(points.size());
    }
};
} // namespace ponca
