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

template<class DataPoint, class Adapter>
class KdTreeKNearestIndexQuery : public KdTreeQuery<DataPoint, Adapter>,
    public KNearestIndexQuery<typename Adapter::IndexType, typename DataPoint::Scalar>
{
public:
    using IndexType       = typename Adapter::IndexType;
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = KNearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType  = KdTreeQuery<DataPoint, Adapter>;

    KdTreeKNearestIndexQuery(const KdTree<DataPoint, Adapter>* kdtree, IndexType k, IndexType index) :
        KdTreeQuery<DataPoint, Adapter>(kdtree), KNearestIndexQuery<IndexType, Scalar>(k, index)
    {
    }

public:
    KdTreeKNearestIterator<IndexType, DataPoint> begin();
    KdTreeKNearestIterator<IndexType, DataPoint> end();

protected:
    void search();
};

#include "./kdTreeKNearestIndexQuery.hpp"
} // namespace ponca
