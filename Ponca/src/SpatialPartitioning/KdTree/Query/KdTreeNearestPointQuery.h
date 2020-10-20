#pragma once

#include <PCA/SpacePartitioning/Query/NearestPointQuery.h>
#include <PCA/SpacePartitioning/KdTree/Query/KdTreeQuery.h>
#include <PCA/SpacePartitioning/KdTree/Iterator/KdTreeNearestPointIterator.h>

namespace pca {

class KdTreeNearestPointQuery : public KdTreeQuery,
                                public NearestPointQuery
{
public:
    KdTreeNearestPointQuery();
    KdTreeNearestPointQuery(const KdTree* kdtree);
    KdTreeNearestPointQuery(const KdTree* kdtree, const Vector3& point);

public:
    KdTreeNearestPointIterator begin();
    KdTreeNearestPointIterator end();

protected:
    void search();
};

} // namespace pca
