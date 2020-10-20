#pragma once


#include "../../query.h"
#include "../query.h"
#include "../kdTree.h"
#include "../iterator.h"

namespace Ponca {

template <typename _VectorType>
struct KdTreeKNearestPointQuery : public KdTreeQuery,
                                 public KNearestPointQuery<_VectorType>
{
    using VectorType = _VectorType;

    KdTreeKNearestPointQuery();
    KdTreeKNearestPointQuery(const KdTree* kdtree, int k);
    KdTreeKNearestPointQuery(const KdTree* kdtree, int k, const VectorType& point);

public:
    KdTreeKNearestPointIterator begin();
    KdTreeKNearestPointIterator end();

protected:
    void search();
};

} // namespace Ponca

#include "./kdTreeKNearestPointQuery.hpp"
