#include <NearestPointQuery.h>

namespace pca {

    NearestPointQuery::NearestPointQuery() :
        PointQuery(),
        NearestQuery()
    {
    }

    NearestPointQuery::NearestPointQuery(const Vector3& point) :
        PointQuery(point),
        NearestQuery()
    {
    }

} // namespace pca