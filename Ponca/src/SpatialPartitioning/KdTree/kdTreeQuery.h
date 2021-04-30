/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../indexSquaredDistance.h"
#include "../../Common/Containers/stack.h"


#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {
template<class DataPoint> class KdTree;

/// \ingroup spatialpartitioning
template <class DataPoint>
class KdTreeQuery
{
public:
    inline KdTreeQuery(const KdTree<DataPoint>* kdtree = nullptr) : m_kdtree( kdtree ), m_stack() {}

protected:
    const KdTree<DataPoint>* m_kdtree { nullptr };
    Stack<IndexSquaredDistance<typename DataPoint::Scalar>, 2 * PCA_KDTREE_MAX_DEPTH> m_stack;
};

} // namespace Ponca