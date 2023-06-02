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

template <typename Traits>
class KdTreeKNearestPointQuery : public KdTreeQuery<Traits>,
    public KNearestPointQuery<typename Traits::IndexType, typename Traits::DataPoint>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = KNearestPointQuery<IndexType, DataPoint>;
    using QueryAccelType = KdTreeQuery<Traits>;

    KdTreeKNearestPointQuery(const KdTreeBase<Traits>* kdtree, IndexType k, const VectorType& point) :
        KdTreeQuery<Traits>(kdtree), KNearestPointQuery<IndexType, DataPoint>(k, point)
    {
    }

public:
    KdTreeKNearestIterator<IndexType, DataPoint> begin();
    KdTreeKNearestIterator<IndexType, DataPoint> end();

protected:
   void search();
};

#include "./kdTreeKNearestPointQuery.hpp"
} // namespace Ponca

