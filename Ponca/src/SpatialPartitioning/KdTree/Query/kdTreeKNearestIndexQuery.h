/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeKNearestIterator.h"

namespace Ponca {

template <typename Traits>
class KdTreeKNearestIndexQuery : public KdTreeQuery<Traits>,
    public KNearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = KNearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeKNearestIterator<IndexType, DataPoint>;

    KdTreeKNearestIndexQuery(const KdTreeBase<Traits>* kdtree, IndexType k, IndexType index) :
        KdTreeQuery<Traits>(kdtree), KNearestIndexQuery<IndexType, Scalar>(k, index)
    {
    }

public:
    Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        this->search();
        return Iterator(QueryType::m_queue.begin());
    }
    Iterator end(){
        return Iterator(QueryType::m_queue.end());
    }

protected:
    void search(){
        KdTreeQuery<Traits>::search_internal(QueryAccelType::m_kdtree->point_data()[QueryType::input()].pos(),
                                             [](IndexType, IndexType){},
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [this](IndexType idx){return QueryType::input() == idx;},
                                             [this](IndexType idx, IndexType, Scalar d){QueryType::m_queue.push({idx, d}); return false;}
        );
    }
};
} // namespace ponca
