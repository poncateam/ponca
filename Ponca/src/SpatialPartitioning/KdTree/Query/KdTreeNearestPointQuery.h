/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeNearestPointIterator.h"
#include "../query.h"


namespace Ponca {

template <class DataPoint>
class KdTreeNearestPointQuery : public NearestPointQuery<DataPoint>, 
								public KdTreeQuery<DataPoint>
{
public:
    using VectorType = typename DataPoint::VectorType;

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
    KdTreeNearestPointIterator begin();
    KdTreeNearestPointIterator end();

protected:
    void search();
};

#include "./KdTreeNearestPointQuery.hpp"
} // namespace ponca
