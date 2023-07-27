/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "knnGraphQuery.h"
#include "../knnGraph.h"

namespace Ponca {

KnnGraphQuery::KnnGraphQuery() :
    m_graph(nullptr),
    m_index(-1)
{
}

KnnGraphQuery::KnnGraphQuery(const KnnGraph* graph, int index) :
    m_graph(graph),
    m_index(index)
{
}

KnnGraphQuery::iterator KnnGraphQuery::begin() const
{
    return m_graph->index_data().begin() + m_index * m_graph->k();
}

KnnGraphQuery::iterator KnnGraphQuery::end() const
{
    return m_graph->index_data().begin() + (m_index+1) * m_graph->k();
}

} // namespace Ponca
