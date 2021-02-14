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

class KdTreeRangeIndexQuery : public KdTreeQuery,
                              public RangeIndexQuery
{
protected:
    friend class KdTreeRangeIndexIterator;

public:
    inline KdTreeRangeIndexQuery();
    inline KdTreeRangeIndexQuery(const KdTree* kdtree);
    inline KdTreeRangeIndexQuery(const KdTree* kdtree, Scalar radius);
    inline KdTreeRangeIndexQuery(const KdTree* kdtree, Scalar radius, int index);

public:
    inline KdTreeRangeIndexIterator begin();
    inline KdTreeRangeIndexIterator end();

protected:
    inline void initialize(KdTreeRangeIndexIterator& iterator);
    inline void advance(KdTreeRangeIndexIterator& iterator);
};

}
#include "./KdTreeRangeIndexQuery.hpp"
