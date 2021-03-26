/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeNearestIndexIterator.h"
#include "../kdTree.h"
#include "../../query.h"

namespace Ponca {

template <class DataPoint>
class KdTreeNearestIndexQuery : public KdTreeQuery<DataPoint>,
								public NearestIndexQuery<typename DataPoint::Scalar>
{
public:
	
    typedef typename DataPoint::Scalar Scalar;

    KdTreeNearestIndexQuery() :
        KdTreeQuery<DataPoint>(), NearestIndexQuery<Scalar>()
    {
    }

    KdTreeNearestIndexQuery(const KdTree<DataPoint>* kdtree) :
        KdTreeQuery<DataPoint>(kdtree), NearestIndexQuery<Scalar>()
    {
    }

    KdTreeNearestIndexQuery(const KdTree<DataPoint>* kdtree, int index) :
        KdTreeQuery<DataPoint>(kdtree), NearestIndexQuery<Scalar>(index)
    {
    }

public:
    KdTreeNearestIndexIterator begin();
    KdTreeNearestIndexIterator end();

protected:
    void search();
};

#include "./KdTreeNearestIndexQuery.hpp"
}
