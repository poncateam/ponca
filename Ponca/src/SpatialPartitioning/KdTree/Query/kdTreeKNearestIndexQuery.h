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
class KdTreeKNearestIndexQuery : public KdTreeQuery<Traits>,
    public KNearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = KNearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;

    KdTreeKNearestIndexQuery(const KdTreeBase<Traits>* kdtree, IndexType k, IndexType index) :
        KdTreeQuery<Traits>(kdtree), KNearestIndexQuery<IndexType, Scalar>(k, index)
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
