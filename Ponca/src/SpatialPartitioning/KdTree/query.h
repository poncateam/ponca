#pragma once

#include "../indexSquaredDistance.h"
#include "./kdTree.h"
#include "../../Common/Containers/stack.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {

class KdTree;

struct KdTreeQuery
{
public:
    inline KdTreeQuery(const KdTree* kdtree = nullptr) : m_kdtree( kdtree ){}

protected:
    const KdTree* m_kdtree { nullptr };
    Stack<IndexSquaredDistance, 2*PCA_KDTREE_MAX_DEPTH> m_stack;
};

//struct KdTreeKNearestQuery;


} // namespace Ponca