/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeNearestIterator.h"

namespace Ponca {

template <typename Traits,
          template <typename> typename IteratorType,
          typename QueryType>
class KdTreeNearestQueryBase : public KdTreeQuery<Traits>, public QueryType
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = IteratorType<typename Traits::IndexType>;

    KdTreeNearestQueryBase(const KdTreeBase<Traits>* kdtree, typename QueryType::InputType input) :
            KdTreeQuery<Traits>(kdtree), QueryType(input){}

public:
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        this->search();
        return Iterator(QueryType::m_nearest);
    }
    inline Iterator end(){
        return Iterator(QueryType::m_nearest + 1);
    }

protected:
    inline void search(){
        KdTreeQuery<Traits>::search_internal(QueryType::getInputPosition(QueryAccelType::m_kdtree->point_data()),
                                             [](IndexType, IndexType){},
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [this](IndexType idx){return QueryType::skipIndexFunctor(idx);},
                                             [this](IndexType idx, IndexType, Scalar d)
                                             {
                                                 QueryType::m_nearest = idx;
                                                 QueryType::m_squared_distance = d;
                                                 return false;
                                             }
        );
    }
};

template <typename Traits>
using KdTreeNearestIndexQuery = KdTreeNearestQueryBase< Traits, KdTreeNearestIterator,
                                NearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>>;
template <typename Traits>
using KdTreeNearestPointQuery = KdTreeNearestQueryBase< Traits, KdTreeNearestIterator,
                                NearestPointQuery<typename Traits::IndexType, typename Traits::DataPoint>>;
} // namespace ponca
