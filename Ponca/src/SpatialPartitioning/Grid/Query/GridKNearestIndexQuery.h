#pragma once

#include <PCA/SpacePartitioning/Query/KNearestIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridQuery.h>
#include <PCA/SpacePartitioning/Grid/Iterator/GridKNearestIndexIterator.h>

namespace Ponca {

class GridKNearestIndexQuery : public GridQuery,
                               public KNearestIndexQuery
{
public:
    GridKNearestIndexQuery();
    GridKNearestIndexQuery(const Grid* grid, int k);
    GridKNearestIndexQuery(const Grid* grid, int k, int index);

public:
    GridKNearestIndexIterator begin();
    GridKNearestIndexIterator end();

protected:
    void search();
};

}   
