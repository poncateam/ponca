#pragma once

#include <PCA/SpacePartitioning/Query/RangePointQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridQuery.h>
#include <PCA/SpacePartitioning/Grid/Iterator/GridRangePointIterator.h>

namespace Ponca {

class GridRangePointQuery : public GridQuery,
                            public RangePointQuery
{
protected:
    friend class GridRangePointIterator;

public:
    GridRangePointQuery();
    GridRangePointQuery(const Grid* grid);
    GridRangePointQuery(const Grid* grid, Scalar radius);
    GridRangePointQuery(const Grid* grid, Scalar radius, const Vector3& point);

public:
    GridRangePointIterator begin();
    GridRangePointIterator end();

protected:
    void initialize(GridRangePointIterator& iterator);
    void advance(GridRangePointIterator& iterator);
};

}   
