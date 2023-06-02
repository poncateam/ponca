/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeRangeIterator.h"

namespace Ponca {

template <typename Traits>
class KdTreeRangeIndexQuery : public KdTreeQuery<Traits>,
    public RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = RangeIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = KdTreeRangeIterator<IndexType, DataPoint, KdTreeRangeIndexQuery>;

protected:
    friend Iterator;

public:

    KdTreeRangeIndexQuery(const KdTreeBase<Traits>* kdtree, Scalar radius, IndexType index) :
        KdTreeQuery<Traits>(kdtree), RangeIndexQuery<IndexType, Scalar>(radius, index)
    {
    }

public:
    inline Iterator begin();
    inline Iterator end();

protected:
    inline void initialize();
    inline void advance(Iterator& iterator);
};

#include "./kdTreeRangeIndexQuery.hpp"
} // namespace ponca
