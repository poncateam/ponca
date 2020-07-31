#include <PCA/SpacePartitioning/Grid/Query/GridKNearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Grid.h>

#include <PCA/Common/Assert.h> //TODO remove this include

namespace pca {

GridKNearestPointQuery::GridKNearestPointQuery() :
    GridQuery(),
    KNearestPointQuery()
{
}

GridKNearestPointQuery::GridKNearestPointQuery(const Grid* grid, int k) :
    GridQuery(grid),
    KNearestPointQuery(k)
{
}

GridKNearestPointQuery::GridKNearestPointQuery(const Grid* grid, int k, const Vector3& point) :
    GridQuery(grid),
    KNearestPointQuery(k, point)
{
}

GridKNearestPointIterator GridKNearestPointQuery::begin()
{
    this->search();
    return GridKNearestPointIterator(m_queue.begin());
}

GridKNearestPointIterator GridKNearestPointQuery::end()
{
    return GridKNearestPointIterator(m_queue.end());
}

void GridKNearestPointQuery::search()
{
    PCA_TODO;
}

} // namespace pca
