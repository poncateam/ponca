/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/kdTreeRangeIndexIterator.h"
//#include "../../query.h"
//#include "../query.h"

namespace Ponca {

//template<class DataPoint> class KdTreeRangeIndexIterator;
template<typename Scalar> class RangeIndexQuery;
template<class DataPoint> class KdTreeQuery;
template<class DataPoint> class KdTree;

template <class DataPoint>
class KdTreeRangeIndexQuery : public KdTreeQuery<DataPoint>, public RangeIndexQuery<typename DataPoint::Scalar>
{
	using VectorType = typename DataPoint::VectorType;
	using Scalar = typename DataPoint::Scalar;

protected:
	friend class KdTreeRangeIndexIterator<DataPoint>;

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
    inline KdTreeRangeIndexIterator<DataPoint> begin();
    inline KdTreeRangeIndexIterator<DataPoint> end();

protected:
    inline void initialize(KdTreeRangeIndexIterator<DataPoint>& iterator);
    inline void advance(KdTreeRangeIndexIterator<DataPoint>& iterator);
};

#include "./kdTreeRangeIndexQuery.hpp"
} // namespace ponca
