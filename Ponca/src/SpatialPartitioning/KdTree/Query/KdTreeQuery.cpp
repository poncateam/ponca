#include <KdTreeQuery.h>

namespace pca {

    KdTreeQuery::KdTreeQuery() :
        m_kdtree(nullptr)
    {
    }

    KdTreeQuery::KdTreeQuery(const KdTree* kdtree) :
        m_kdtree(kdtree)
    {
    }

} // namespace pca