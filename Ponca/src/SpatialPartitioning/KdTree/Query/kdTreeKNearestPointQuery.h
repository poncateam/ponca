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
class KdTreeKNearestPointQuery : public KdTreeQuery<Traits>,
    public KNearestPointQuery<typename Traits::IndexType, typename Traits::DataPoint>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = KNearestPointQuery<IndexType, DataPoint>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeKNearestIterator<IndexType, DataPoint>;

    inline KdTreeKNearestPointQuery(const KdTreeBase<Traits>* kdtree, IndexType k, const VectorType& point) :
        KdTreeQuery<Traits>(kdtree), KNearestPointQuery<IndexType, DataPoint>(k, point)
    {
    }

public:
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        this->search();
        return Iterator(QueryType::m_queue.begin());
    }
    inline Iterator end(){
        return Iterator(QueryType::m_queue.end());
    }

protected:
   inline void search(){
       KdTreeQuery<Traits>::search_internal(QueryType::input(),
                                            [](IndexType, IndexType){},
                                            [this](){return QueryType::descentDistanceThreshold();},
                                            [](IndexType){return false;},
                                            [this](IndexType idx, IndexType, Scalar d){QueryType::m_queue.push({idx, d}); return false;}
       );
   }
};
} // namespace Ponca
