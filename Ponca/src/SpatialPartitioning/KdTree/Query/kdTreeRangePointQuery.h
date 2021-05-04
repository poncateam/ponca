/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeRangePointIterator.h"

namespace Ponca {

template <class DataPoint>
class KdTreeRangePointQuery : public KdTreeQuery<DataPoint>, public RangePointQuery<DataPoint>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = RangePointQuery<DataPoint>;
    using QueryAccelType  = KdTreeQuery<DataPoint>;

protected:
    friend class KdTreeRangePointIterator<DataPoint>;

public:
    inline KdTreeRangePointQuery() :
        KdTreeQuery<DataPoint>(), RangePointQuery<DataPoint>()
    {
    }

    inline KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>()
    {
    }

    inline KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>(radius)
    {
    }

    inline KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius, const VectorType& point) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>(radius, point)
    {
    }

public:
    inline KdTreeRangePointIterator<DataPoint> begin();
	inline KdTreeRangePointIterator<DataPoint> end();

protected:
	inline void initialize(KdTreeRangePointIterator<DataPoint>& iterator);
	inline void advance(KdTreeRangePointIterator<DataPoint>& iterator);
};

#include "./kdTreeRangePointQuery.hpp"
} // namespace ponca
