/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeKNearestPointIterator.h"

#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"

namespace Ponca {

template <class DataPoint>
class KdTreeKNearestPointQuery : public KNearestPointQuery<typename DataPoint::VectorType>, public KdTreeQuery<DataPoint>
{
public:
    using VectorType = DataPoint::VectorType;
    using Scalar = DataPoint::Scalar;

    KdTreeKNearestPointQuery() :
        KdTreeQuery<DataPoint>(), KNearestPointQuery<DataPoint>()
    {
        cout << "Test" << endl;
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
    KdTreeKNearestPointIterator begin();
    KdTreeKNearestPointIterator end();

protected:
   void search();
};

#include "./KdTreeKNearestPointQuery.hpp"
} // namespace Ponca

