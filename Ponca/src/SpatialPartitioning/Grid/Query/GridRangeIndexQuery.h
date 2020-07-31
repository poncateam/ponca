#pragma once

#include <PCA/SpacePartitioning/Query/RangeIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridQuery.h>
#include <PCA/SpacePartitioning/Grid/Iterator/GridRangeIndexIterator.h>

namespace pca {

class GridRangeIndexQuery : public GridQuery,
                            public RangeIndexQuery
{
protected:
    friend class GridRangeIndexIterator;

public:
    GridRangeIndexQuery();
    GridRangeIndexQuery(const Grid* grid);
    GridRangeIndexQuery(const Grid* grid, Scalar radius);
    GridRangeIndexQuery(const Grid* grid, Scalar radius, int index);

public:
    GridRangeIndexIterator begin();
    GridRangeIndexIterator end();

protected:
    void initialize(GridRangeIndexIterator& iterator);
    void advance(GridRangeIndexIterator& iterator);
};

} // namespace pca
