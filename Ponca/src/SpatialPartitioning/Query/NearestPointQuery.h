#pragma once

#include <Ponca/src/SpatialPartitioning/Query/PointQuery.h>
#include <Ponca/src/SpatialPartitioning/Query/NearestQuery.h>

namespace pca {

    class NearestPointQuery : public PointQuery,
        public NearestQuery
    {
    public:
        NearestPointQuery();
        NearestPointQuery(const Vector3& point);
    };

} // namespace pca