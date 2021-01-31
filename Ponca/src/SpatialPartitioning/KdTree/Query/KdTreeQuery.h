#pragma once

#include "../../indexSquaredDistance.h"
#include "../../../Common/Containers/stack.h"

#define PDPC_KDTREE_MAX_DEPTH 32

namespace pca {

    class KdTree;

    class KdTreeQuery
    {
    public:
        KdTreeQuery();
        KdTreeQuery(const KdTree* kdtree);

    protected:
        const KdTree* m_kdtree;
        Stack<IndexSquaredDistance, 2 * PDPC_KDTREE_MAX_DEPTH> m_stack;
    };

} // namespace pca