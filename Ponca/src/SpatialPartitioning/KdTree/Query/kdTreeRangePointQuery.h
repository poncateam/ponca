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
class KdTreeRangePointQuery : public KdTreeQuery<Traits>,
    public RangePointQuery<typename Traits::IndexType, typename Traits::DataPoint>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = RangePointQuery<IndexType, DataPoint>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeRangeIterator<IndexType, DataPoint, KdTreeRangePointQuery>;

protected:
    friend Iterator;

public:

    inline KdTreeRangePointQuery(const KdTreeBase<Traits>* kdtree, Scalar radius, const VectorType& point) :
        KdTreeQuery<Traits>(kdtree), RangePointQuery<IndexType, DataPoint>(radius, point)
    {
    }

public:
    inline Iterator begin();
    inline Iterator end();

protected:
    inline void advance(Iterator& iterator);
};

#include "./kdTreeRangePointQuery.hpp"
} // namespace ponca
