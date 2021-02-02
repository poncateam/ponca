#pragma once

#include <PCA/SpacePartitioning/Query/NearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridQuery.h>
#include <PCA/SpacePartitioning/Grid/Iterator/GridNearestPointIterator.h>

namespace Ponca {

class GridNearestPointQuery : public GridQuery,
                              public NearestPointQuery
{
public:
    GridNearestPointQuery();
    GridNearestPointQuery(const Grid* grid);
    GridNearestPointQuery(const Grid* grid, const Vector3& point);

public:
    GridNearestPointIterator begin();
    GridNearestPointIterator end();

protected:
    void search();
};

}   
