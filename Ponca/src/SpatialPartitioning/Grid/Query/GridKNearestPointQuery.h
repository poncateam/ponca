#pragma once

#include <PCA/SpacePartitioning/Query/KNearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridQuery.h>
#include <PCA/SpacePartitioning/Grid/Iterator/GridKNearestPointIterator.h>

namespace pca {

class GridKNearestPointQuery : public GridQuery,
                               public KNearestPointQuery
{
public:
    GridKNearestPointQuery();
    GridKNearestPointQuery(const Grid* grid, int k);
    GridKNearestPointQuery(const Grid* grid, int k, const Vector3& point);

public:
    GridKNearestPointIterator begin();
    GridKNearestPointIterator end();

protected:
    void search();
};

} // namespace pca
