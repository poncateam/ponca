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
template <class DataPoint, class Compatibility> class KdTree;

template <class DataPoint, class Compatibility>
class KdTreeQuery
{
public:
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;

    explicit inline KdTreeQuery(const KdTree<DataPoint, Compatibility>* kdtree) : m_kdtree( kdtree ), m_stack() {}

protected:
    /// \brief Init stack for a new search
    inline void reset() {
        m_stack.clear();
        m_stack.push({0,0});
    }

    const KdTree<DataPoint, Compatibility>* m_kdtree { nullptr };
    Stack<IndexSquaredDistance<typename DataPoint::Scalar>, 2 * PCA_KDTREE_MAX_DEPTH> m_stack;
};

} // namespace Ponca