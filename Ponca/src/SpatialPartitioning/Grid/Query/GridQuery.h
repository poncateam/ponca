#pragma once

namespace Ponca {

class Grid;

class GridQuery
{
public:
    GridQuery();
    GridQuery(const Grid* grid);

protected:
    const Grid* m_grid;
};

}   
