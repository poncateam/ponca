/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeRangePointIterator.h"
//
//#include "../query.h"
//#include "../../query.h"

namespace Ponca {

template<class DataPoint> class RangePointQuery;
template<class DataPoint> class KdTreeQuery;
template<class DataPoint> class KdTree;
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
        KdTreeQuery<DataPoint>(kdtree), RangePointQuery<VectorType>(radius, point)
    {
    }

public:
    inline KdTreeRangePointIterator<DataPoint> begin();
	inline KdTreeRangePointIterator<DataPoint> end();

protected:
	inline void initialize(KdTreeRangePointIterator<DataPoint>& iterator);
	inline void advance(KdTreeRangePointIterator<DataPoint>& iterator);
};

#include "./KdTreeRangePointQuery.hpp"
} // namespace ponca
