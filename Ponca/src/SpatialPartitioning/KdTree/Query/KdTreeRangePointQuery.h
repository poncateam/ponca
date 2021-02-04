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

template <typename _VectorType>
class KdTreeRangePointQuery : public RangePointQuery<_VectorType>, public KdTreeQuery
{
using VectorType = typename RangePointQuery<_VectorType>::VectorType;

protected:
    friend class KdTreeRangePointIterator<VectorType>;

public:
    KdTreeRangePointQuery();
    KdTreeRangePointQuery(const KdTree* kdtree);
    KdTreeRangePointQuery(const KdTree* kdtree, Scalar radius);
    KdTreeRangePointQuery(const KdTree* kdtree, Scalar radius, const VectorType& point);

public:
    KdTreeRangePointIterator<VectorType> begin();
    KdTreeRangePointIterator<VectorType> end();

protected:
    void initialize(KdTreeRangePointIterator<VectorType>& iterator);
    void advance(KdTreeRangePointIterator<VectorType>& iterator);
};

} // namespace ponca
#include "./KdTreeRangePointQuery.hpp"
