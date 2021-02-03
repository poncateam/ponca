/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/KdTreeNearestPointIterator.h"

#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"

namespace Ponca {

template <typename _VectorType>
class KdTreeNearestPointQuery : public NearestPointQuery<_VectorType>, public KdTreeQuery
{
using VectorType = typename NearestPointQuery<_VectorType>::VectorType;
public:
    KdTreeNearestPointQuery();
    KdTreeNearestPointQuery(const KdTree* kdtree);
    KdTreeNearestPointQuery(const KdTree* kdtree, const VectorType& point);

public:
    KdTreeNearestPointIterator begin();
    KdTreeNearestPointIterator end();

protected:
    void search();
};

}   
#include "./KdTreeNearestPointQuery.hpp"
