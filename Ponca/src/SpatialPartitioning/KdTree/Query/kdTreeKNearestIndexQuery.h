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

template<class DataPoint, class Compatibility>
class KdTreeKNearestIndexQuery : public KdTreeQuery<DataPoint, Compatibility>,
    public KNearestIndexQuery<typename DataPoint::Scalar>
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = KNearestIndexQuery<typename DataPoint::Scalar>;
    using QueryAccelType  = KdTreeQuery<DataPoint, Compatibility>;

    KdTreeKNearestIndexQuery(const KdTree<DataPoint, Compatibility>* kdtree, int k, int index) :
        KdTreeQuery<DataPoint, Compatibility>(kdtree), KNearestIndexQuery<Scalar>(k, index)
    {
    }

public:
    KdTreeKNearestIterator<DataPoint> begin();
    KdTreeKNearestIterator<DataPoint> end();

protected:
    void search();
};

#include "./kdTreeKNearestIndexQuery.hpp"
} // namespace ponca
