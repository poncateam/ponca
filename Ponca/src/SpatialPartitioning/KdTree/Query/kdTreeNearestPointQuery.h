/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeNearestIterator.h"

namespace Ponca {

template <class DataPoint>
class KdTreeNearestPointQuery : public NearestPointQuery<DataPoint>, public KdTreeQuery<DataPoint>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = NearestPointQuery<DataPoint>;
    using QueryAccelType  = KdTreeQuery<DataPoint>;

    KdTreeNearestPointQuery() :
        KdTreeQuery<DataPoint>(), NearestPointQuery<DataPoint>()
    {
    }

    KdTreeNearestPointQuery(const KdTree<DataPoint>* kdtree) :
        KdTreeQuery<DataPoint>(kdtree), NearestPointQuery<DataPoint>()
    {
    }

    KdTreeNearestPointQuery(const KdTree<DataPoint>* kdtree, const VectorType& point) :
        KdTreeQuery<DataPoint>(kdtree), NearestPointQuery<DataPoint>(point)
    {
    }

public:
    KdTreeNearestIterator begin();
    KdTreeNearestIterator end();

protected:
    void search();
};

#include "./kdTreeNearestPointQuery.hpp"
} // namespace ponca
