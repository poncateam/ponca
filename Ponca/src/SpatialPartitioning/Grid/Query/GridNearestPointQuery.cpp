#include <PCA/SpacePartitioning/Grid/Query/GridNearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Grid.h>

#include <PCA/Common/Assert.h> //TODO remove this include

namespace Ponca {

GridNearestPointQuery::GridNearestPointQuery() :
    GridQuery(),
    NearestPointQuery()
{
}

GridNearestPointQuery::GridNearestPointQuery(const Grid* grid) :
    GridQuery(grid),
    NearestPointQuery()
{
}

GridNearestPointQuery::GridNearestPointQuery(const Grid* grid, const Vector3& point) :
    GridQuery(grid),
    NearestPointQuery(point)
{
}

GridNearestPointIterator GridNearestPointQuery::begin()
{
    this->search();
    return GridNearestPointIterator(m_nearest);
}

GridNearestPointIterator GridNearestPointQuery::end()
{
    return GridNearestPointIterator(m_nearest+1);
}

void GridNearestPointQuery::search()
{
    PCA_TODO;
}

} // namespace pca
