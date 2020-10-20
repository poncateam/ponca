#pragma once

#include <PCA/SpacePartitioning/Query/NearestIndexQuery.h>
#include <PCA/SpacePartitioning/KdTree/Query/KdTreeQuery.h>
#include <PCA/SpacePartitioning/KdTree/Iterator/KdTreeNearestIndexIterator.h>

namespace pca {

class KdTreeNearestIndexQuery : public KdTreeQuery,
                                public NearestIndexQuery
{
public:
    KdTreeNearestIndexQuery();
    KdTreeNearestIndexQuery(const KdTree* kdtree);
    KdTreeNearestIndexQuery(const KdTree* kdtree, int index);

public:
    KdTreeNearestIndexIterator begin();
    KdTreeNearestIndexIterator end();

protected:
    void search();
};

} // namespace pca
