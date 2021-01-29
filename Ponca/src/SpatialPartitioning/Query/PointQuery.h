#pragma once

#include <Ponca/src/Common/Defines.h>

namespace pca {

    class PointQuery
    {
    public:
        PointQuery();
        PointQuery(const Vector3& point);

        const Vector3& point() const;
        void set_point(const Vector3& point);

    protected:
        Vector3 m_point;
    };

} // namespace pca