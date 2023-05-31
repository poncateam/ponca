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

template <class DataPoint, class Adapter>
class KdTreeNearestIndexQuery : public KdTreeQuery<DataPoint, Adapter>,
    public NearestIndexQuery<typename Adapter::IndexType, typename DataPoint::Scalar>
{
public:
    using IndexType       = typename Adapter::IndexType;
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = NearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType  = KdTreeQuery<DataPoint, Adapter>;

    KdTreeNearestIndexQuery(const KdTree<DataPoint, Adapter>* kdtree, IndexType index) :
        KdTreeQuery<DataPoint, Adapter>(kdtree), NearestIndexQuery<IndexType, Scalar>(index)
    {
    }

public:
    KdTreeNearestIterator<IndexType> begin();
    KdTreeNearestIterator<IndexType> end();

protected:
    void search();
};

#include "./kdTreeNearestIndexQuery.hpp"
} // namespace ponca
