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
	using VectorType = typename DataPoint::VectorType;
	using Scalar = typename DataPoint::Scalar;

protected:
    friend class KdTreeRangePointIterator<DataPoint>;

public:
    KdTreeRangePointQuery() :
        KdTreeQuery<DataPoint>(), RangePointQuery<DataPoint>()
    {
    }

    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>()
    {
    }

    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius) :
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<DataPoint>(radius)
    {
    }

    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius, const VectorType& point) :
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
