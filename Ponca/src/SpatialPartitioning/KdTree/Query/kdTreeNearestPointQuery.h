/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeNearestIterator.h"

namespace Ponca {

template <typename Traits>
class KdTreeNearestPointQuery : public KdTreeQuery<Traits>,
    public NearestPointQuery<typename Traits::IndexType, typename Traits::DataPoint>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = NearestPointQuery<IndexType, DataPoint>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeNearestIterator<IndexType>;

    inline KdTreeNearestPointQuery(const KdTreeBase<Traits>* kdtree, const VectorType& point) :
        KdTreeQuery<Traits>(kdtree), NearestPointQuery<IndexType, DataPoint>(point)
    {
    }

public:
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        this->search();
        return Iterator(QueryType::m_nearest);
    }
    Iterator end(){
        return Iterator(QueryType::m_nearest + 1);
    }

protected:
    inline void search(){
        KdTreeQuery<Traits>::search_internal(QueryType::input(),
                                             [](IndexType, IndexType){},
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [](IndexType){return false;},
                                             [this](IndexType idx, IndexType, Scalar d)
                                             {
                                                 QueryType::m_nearest = idx;
                                                 QueryType::m_squared_distance = d;
                                                 return false;
                                             }
        );
    }
};
} // namespace ponca
