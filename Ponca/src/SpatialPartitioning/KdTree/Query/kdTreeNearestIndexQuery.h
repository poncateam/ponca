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

template <typename Traits>
class KdTreeNearestIndexQuery : public KdTreeQuery<Traits>,
    public NearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryType      = NearestIndexQuery<IndexType, typename DataPoint::Scalar>;
    using QueryAccelType = KdTreeQuery<Traits>;

    KdTreeNearestIndexQuery(const KdTreeBase<Traits>* kdtree, IndexType index) :
        KdTreeQuery<Traits>(kdtree), NearestIndexQuery<IndexType, Scalar>(index)
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
