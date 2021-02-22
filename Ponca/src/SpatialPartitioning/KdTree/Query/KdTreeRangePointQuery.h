/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeRangePointIterator.h"

#include "../../query.h"
#include "../query.h"

namespace Ponca {

template <class DataPoint>
class KdTreeRangePointQuery : public RangePointQuery<DataPoint>, public KdTreeQuery<DataPoint>
{
	using VectorType = typename DataPoint::VectorType;
	using Scalar = typename DataPoint::Scalar;

protected:
    friend class KdTreeRangePointIterator<VectorType>;

public:
    KdTreeRangePointQuery();
    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree);
    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius);
    KdTreeRangePointQuery(const KdTree<DataPoint>* kdtree, Scalar radius, const VectorType& point);

public:
    KdTreeRangePointIterator<VectorType> begin();
    KdTreeRangePointIterator<VectorType> end();

protected:
    void initialize(KdTreeRangePointIterator<VectorType>& iterator);
    void advance(KdTreeRangePointIterator<VectorType>& iterator);
};

} // namespace ponca
#include "./KdTreeRangePointQuery.hpp"
