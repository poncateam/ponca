/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeRangeIterator.h"

namespace Ponca {

/*!
 * \brief Extension of the Query class that allows to read the result of a range neighbors search on the KdTree.
 *
 *  Ouput result of a `KdTree::rangeNeighbors` query request.
 *
 *  \see KdTreeBase
 */
template <typename Traits,
        template <typename,typename,typename> typename IteratorType,
        typename QueryType>
class KdTreeRangeQueryBase : public KdTreeQuery<Traits>, public QueryType
                              //public RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = IteratorType<IndexType, DataPoint, KdTreeRangeQueryBase>;
    using Self           = KdTreeRangeQueryBase<Traits, IteratorType, QueryType>;

protected:
    friend Iterator;

public:
    KdTreeRangeQueryBase(const KdTreeBase<Traits>* kdtree, Scalar radius, typename QueryType::InputType input) :
            KdTreeQuery<Traits>(kdtree), QueryType(radius, input){}

    /// \brief Call the range neighbors query with new input and radius parameters.
    inline Self& operator()(typename QueryType::InputType input, Scalar radius)
    {
        return QueryType::template operator()<Self>(input, radius);
    }

    /// \brief Call the range neighbors query with new input parameter.
    inline Self& operator()(typename QueryType::InputType input)
    {
        return QueryType::template operator()<Self>(input);
    }

    /// \brief Returns an iterator to the beginning of the Range Query.
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        Iterator it(this);
        this->advance(it);
        return it;
    }

    /// \brief Returns an iterator to the end of the Range Query.
    inline Iterator end(){
        return Iterator(this, QueryAccelType::m_kdtree->pointCount());
    }

protected:
    inline void advance(Iterator& it){
        const auto& points  = QueryAccelType::m_kdtree->points();
        const auto& indices = QueryAccelType::m_kdtree->samples();
        const auto& point   = QueryType::getInputPosition(points);

        if (points.empty() || indices.empty())
        {
            it = end();
            return;
        }

        auto descentDistanceThreshold = [this](){return QueryType::descentDistanceThreshold();};
        auto skipFunctor              = [this](IndexType idx){return QueryType::skipIndexFunctor(idx);};
        auto processNeighborFunctor   = [&it](IndexType idx, IndexType i, Scalar)
        {
            it.m_index = idx;
            it.m_start = i+1;
            return true;
        };

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

template <typename Traits>
using KdTreeRangeIndexQuery = KdTreeRangeQueryBase< Traits, KdTreeRangeIterator,
        RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>>;
template <typename Traits>
using KdTreeRangePointQuery = KdTreeRangeQueryBase< Traits, KdTreeRangeIterator,
        RangePointQuery<typename Traits::IndexType, typename Traits::DataPoint>>;
} // namespace ponca
