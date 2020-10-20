#include <PCA/SpacePartitioning/Grid/Query/GridKNearestIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Grid.h>

#include <PCA/Common/Assert.h> //TODO remove this include

namespace pca {

GridKNearestIndexQuery::GridKNearestIndexQuery() :
    GridQuery(),
    KNearestIndexQuery()
{
}

GridKNearestIndexQuery::GridKNearestIndexQuery(const Grid* grid, int k) :
    GridQuery(grid),
    KNearestIndexQuery(k)
{
}

GridKNearestIndexQuery::GridKNearestIndexQuery(const Grid* grid, int k, int index) :
    GridQuery(grid),
    KNearestIndexQuery(k, index)
{
}

GridKNearestIndexIterator GridKNearestIndexQuery::begin()
{
    this->search();
    return GridKNearestIndexIterator(m_queue.begin());
}

GridKNearestIndexIterator GridKNearestIndexQuery::end()
{
    return GridKNearestIndexIterator(m_queue.end());
}

void GridKNearestIndexQuery::search()
{
    PCA_TODO;
}

} // namespace pca
