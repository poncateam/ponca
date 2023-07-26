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
class KdTreeNearestIndexQuery : public KdTreeQuery<Traits>,
    public NearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = NearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeNearestIterator<IndexType>;

    KdTreeNearestIndexQuery(const KdTreeBase<Traits>* kdtree, IndexType index) :
        KdTreeQuery<Traits>(kdtree), NearestIndexQuery<IndexType, Scalar>(index)
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
    void search(){
        KdTreeQuery<Traits>::search_internal(QueryAccelType::m_kdtree->point_data()[QueryType::input()].pos(),
                                             [](IndexType, IndexType){},
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [this](IndexType idx){return QueryType::input() == idx;},
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
