#include <PCA/SpacePartitioning/Grid/Query/GridNearestIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Grid.h>

#include <PCA/Common/Assert.h> //TODO remove this include

namespace Ponca {

GridNearestIndexQuery::GridNearestIndexQuery() :
    GridQuery(),
    NearestIndexQuery()
{
}

GridNearestIndexQuery::GridNearestIndexQuery(const Grid* grid) :
    GridQuery(grid),
    NearestIndexQuery()
{
}

GridNearestIndexQuery::GridNearestIndexQuery(const Grid* grid, int index) :
    GridQuery(grid),
    NearestIndexQuery(index)
{
}

GridNearestIndexIterator GridNearestIndexQuery::begin()
{
    this->search();
    return GridNearestIndexIterator(m_nearest);
}

GridNearestIndexIterator GridNearestIndexQuery::end()
{
    return GridNearestIndexIterator(m_nearest+1);
}

void GridNearestIndexQuery::search()
{
    PCA_TODO;
}

}   
