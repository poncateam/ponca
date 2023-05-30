/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeKNearestIterator.h"

namespace Ponca {

template <class DataPoint, class Adapter>
class KdTreeKNearestPointQuery : public KNearestPointQuery<DataPoint>, public KdTreeQuery<DataPoint, Adapter>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = KNearestPointQuery<DataPoint>;
    using QueryAccelType  = KdTreeQuery<DataPoint, Adapter>;

    KdTreeKNearestPointQuery(const KdTree<DataPoint, Adapter>* kdtree, int k, const VectorType& point) :
        KdTreeQuery<DataPoint, Adapter>(kdtree), KNearestPointQuery<DataPoint>(k, point)
    {
    }

public:
    KdTreeKNearestIterator<DataPoint> begin();
    KdTreeKNearestIterator<DataPoint> end();

protected:
   void search();
};

#include "./kdTreeKNearestPointQuery.hpp"
} // namespace Ponca

