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
class KdTreeRangeIndexQuery : public KdTreeQuery<DataPoint>, public RangeIndexQuery<typename DataPoint::Scalar>
{
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = RangeIndexQuery<typename DataPoint::Scalar>;
    using QueryAccelType  = KdTreeQuery<DataPoint>;
    using Iterator        = KdTreeRangeIterator<DataPoint, KdTreeRangeIndexQuery>;

protected:
	friend Iterator;

public:
    KdTreeRangeIndexQuery() :
        KdTreeQuery<DataPoint>(), RangeIndexQuery<Scalar>()
    {
    }

    KdTreeRangeIndexQuery(const KdTree<DataPoint>* kdtree) :
        KdTreeQuery<DataPoint>(kdtree), RangeIndexQuery<Scalar>()
    {
    }

    KdTreeRangeIndexQuery(const KdTree<DataPoint>* kdtree, Scalar radius) :
        KdTreeQuery<DataPoint>(kdtree), RangeIndexQuery<Scalar>(radius)
    {
    }

    KdTreeRangeIndexQuery(const KdTree<DataPoint>* kdtree, Scalar radius, int index) :
        KdTreeQuery<DataPoint>(kdtree), RangeIndexQuery<Scalar>(radius, index)
    {
    }

public:
    inline Iterator begin();
    inline Iterator end();

protected:
    inline void initialize(Iterator& iterator);
    inline void advance(Iterator& iterator);
};

#include "./kdTreeRangeIndexQuery.hpp"
} // namespace ponca
