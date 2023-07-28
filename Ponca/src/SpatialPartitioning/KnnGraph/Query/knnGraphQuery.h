/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Iterator/knnGraphRangeIterator.h"
#include <vector>

namespace Ponca {

template <typename Traits>class KnnGraphBase; // Need forward declaration to avoid mutual inclusion

/// \todo Inherit from Base queries
template <typename Traits>class KnnGraphQuery
{
public:
    using Iterator = typename Traits::IndexContainer::const_iterator;

public:
    inline KnnGraphQuery(const KnnGraphBase<Traits>* graph, int index)
        : m_graph(graph), m_index(index){}

    inline Iterator begin() const{
        return m_graph->index_data().begin() + m_index * m_graph->k();
    }
    inline Iterator end() const{
        return m_graph->index_data().begin() + (m_index+1) * m_graph->k();
    }

protected:
    const KnnGraphBase<Traits>* m_graph {nullptr};
    int m_index {-1};
};

} // namespace Ponca
