/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../indexSquaredDistance.h"
#include "../../Common/Containers/stack.h"

namespace Ponca {
template <typename Traits> class KdTreeBase;

template <typename Traits>
class KdTreeQuery
{
public:
    using DataPoint  = typename Traits::DataPoint;
    using IndexType  = typename Traits::IndexType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;

    explicit inline KdTreeQuery(const KdTreeBase<Traits>* kdtree) : m_kdtree( kdtree ), m_stack() {}

protected:
    /// \brief Init stack for a new search
    inline void reset() {
        m_stack.clear();
        m_stack.push({0,0});
    }

    const KdTreeBase<Traits>* m_kdtree { nullptr };
    Stack<IndexSquaredDistance<IndexType, Scalar>, 2 * Traits::MAX_DEPTH> m_stack;
};
} // namespace Ponca
