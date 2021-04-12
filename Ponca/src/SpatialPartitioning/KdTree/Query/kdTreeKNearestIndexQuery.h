/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/kdTreeKNearestIndexIterator.h"
#include "../query.h"

namespace Ponca {

template<class DataPoint>
class KdTreeKNearestIndexQuery : public KdTreeQuery<DataPoint>,
    public KNearestIndexQuery<typename DataPoint::Scalar>
{
public:
    typedef typename DataPoint::Scalar Scalar;

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
