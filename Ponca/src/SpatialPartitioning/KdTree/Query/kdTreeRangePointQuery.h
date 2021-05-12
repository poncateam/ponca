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

template <class DataPoint>
class KdTreeRangePointQuery : public KdTreeQuery<DataPoint>, public RangePointQuery<DataPoint>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = RangePointQuery<DataPoint>;
    using QueryAccelType  = KdTreeQuery<DataPoint>;
    using Iterator        = KdTreeRangeIterator<DataPoint, KdTreeRangePointQuery>;


protected:
    friend Iterator;

public:

    inline KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius, const VectorType& point) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>(radius, point)
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
