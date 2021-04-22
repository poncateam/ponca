/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeKNearestPointIterator.h"

namespace Ponca {

template <class DataPoint>
class KdTreeKNearestPointQuery : public KNearestPointQuery<DataPoint>, public KdTreeQuery<DataPoint>
{
public:
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::Scalar Scalar;

    KdTreeKNearestPointQuery() :
        KdTreeQuery<DataPoint>(), KNearestPointQuery<DataPoint>()
    {
    }

    KdTreeKNearestPointQuery(const KdTree<DataPoint>* kdtree, int k) :
        KdTreeQuery<DataPoint>(kdtree), KNearestPointQuery<DataPoint>(k)
    {
    }

    KdTreeKNearestPointQuery(const KdTree<DataPoint>* kdtree, int k, const VectorType& point) :
        KdTreeQuery<DataPoint>(kdtree), KNearestPointQuery<DataPoint>(k, point)
    {
    }

public:
    KdTreeKNearestPointIterator<DataPoint> begin();
    KdTreeKNearestPointIterator<DataPoint> end();

protected:
   void search();
};

#include "./kdTreeKNearestPointQuery.hpp"
} // namespace Ponca

