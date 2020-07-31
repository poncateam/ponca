#pragma once

#include "../iterator.h"
#include "./query.h"

namespace Ponca {

using KdTreeNearestIndexIterator = NearestIterator;
using KdTreeNearestPointIterator = NearestIterator;

using KdTreeKNearestIndexIterator = KNearestIterator;
using KdTreeKNearestPointIterator = KNearestIterator;
//struct KdTreeRangeIndexIterator : public RangeIterator<kdTreeRangeIndexQuery> {};
//struct KdTreeRangePointIterator : public RangeIterator<kdTreeRangePointQuery> {};

} // namespace Ponca