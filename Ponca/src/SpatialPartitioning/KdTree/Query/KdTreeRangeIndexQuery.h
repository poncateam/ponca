/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeRangeIndexIterator.h"
#include "../../query.h"
#include "../query.h"

namespace Ponca {

//
//:	public KdTreeQuery<DataPoint>,
//public RangeIndexQuery<DataPoint::Scalar>
template <class DataPoint>
class KdTreeRangeIndexQuery 
{
	/*using VectorType = typename DataPoint::VectorType;
	using Scalar = typename DataPoint::Scalar;*/
//
protected:
	//template<DataPoint> friend class KdTreeRangeIndexIterator;

public:
  /*  inline KdTreeRangeIndexQuery();
    inline KdTreeRangeIndexQuery(const KdTree* kdtree);
    inline KdTreeRangeIndexQuery(const KdTree* kdtree, Scalar radius);
    inline KdTreeRangeIndexQuery(const KdTree* kdtree, Scalar radius, int index);*/
//
public:
    inline KdTreeRangeIndexIterator<DataPoint> begin();
    inline KdTreeRangeIndexIterator<DataPoint> end();
//
//protected:
//    inline void initialize(KdTreeRangeIndexIterator& iterator);
//    inline void advance(KdTreeRangeIndexIterator& iterator);
};

}
#include "./KdTreeRangeIndexQuery.hpp"
