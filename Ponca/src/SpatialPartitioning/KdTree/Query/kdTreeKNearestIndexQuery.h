/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeKNearestIndexIterator.h"

namespace Ponca {

template<class DataPoint>
class KdTreeKNearestIndexQuery : public KdTreeQuery<DataPoint>,
    public KNearestIndexQuery<typename DataPoint::Scalar>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = KNearestIndexQuery<typename DataPoint::Scalar>;
    using QueryAccelType  = KdTreeQuery<DataPoint>;

    KdTreeKNearestIndexQuery() :
        KdTreeQuery<DataPoint>(), KNearestIndexQuery<Scalar>()
    {
    }

    KdTreeKNearestIndexQuery(const KdTree<DataPoint>* kdtree, int k) :
        KdTreeQuery<DataPoint>(kdtree), KNearestIndexQuery<Scalar>(k)
    {
    }

    KdTreeKNearestIndexQuery(const KdTree<DataPoint>* kdtree, int k, int index) :
        KdTreeQuery<DataPoint>(kdtree), KNearestIndexQuery<Scalar>(k, index)
    {
    }

public:
    KdTreeKNearestIndexIterator<DataPoint> begin();
    KdTreeKNearestIndexIterator<DataPoint> end();

protected:
    void search();
};

#include "./kdTreeKNearestIndexQuery.hpp"
} // namespace ponca
