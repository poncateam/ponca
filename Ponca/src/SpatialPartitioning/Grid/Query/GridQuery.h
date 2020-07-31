#pragma once

namespace pca {

class Grid;

class GridQuery
{
public:
    GridQuery();
    GridQuery(const Grid* grid);

protected:
    const Grid* m_grid;
};

} // namespace pca
