#pragma once

#include "../../query.h"
#include "../query.h"
#include "../iterator.h"

namespace Ponca {

class KdTreeKNearestIndexQuery : public KdTreeQuery,
                                 public KNearestIndexQuery
{
public:
    KdTreeKNearestIndexQuery();
    KdTreeKNearestIndexQuery(const KdTree* kdtree, int k);
    KdTreeKNearestIndexQuery(const KdTree* kdtree, int k, int index);

public:
    KdTreeKNearestIndexIterator begin();
    KdTreeKNearestIndexIterator end();

protected:
    void search();
};

} // namespace pca


#include "./kdTreeKNearestIndexQuery.hpp"