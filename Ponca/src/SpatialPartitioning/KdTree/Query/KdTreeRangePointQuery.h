/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <PCA/SpacePartitioning/Query/RangePointQuery.h>
#include <PCA/SpacePartitioning/KdTree/Query/KdTreeQuery.h>
#include <PCA/SpacePartitioning/KdTree/Iterator/KdTreeRangePointIterator.h>

namespace Ponca {

class KdTreeRangePointQuery : public KdTreeQuery,
                              public RangePointQuery
{
protected:
    friend class KdTreeRangePointIterator;

public:
    KdTreeRangePointQuery();
    KdTreeRangePointQuery(const KdTree* kdtree);
    KdTreeRangePointQuery(const KdTree* kdtree, Scalar radius);
    KdTreeRangePointQuery(const KdTree* kdtree, Scalar radius, const Vector3& point);

public:
    KdTreeRangePointIterator begin();
    KdTreeRangePointIterator end();

protected:
    void initialize(KdTreeRangePointIterator& iterator);
    void advance(KdTreeRangePointIterator& iterator);
};

}   
