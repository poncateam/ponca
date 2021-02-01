#include <PCA/SpacePartitioning/KNNGraph/Query/KNNGraphQuery.h>
#include <PCA/SpacePartitioning/KNNGraph/KNNGraph.h>

namespace Ponca {

KNNGraphQuery::KNNGraphQuery() :
    m_graph(nullptr),
    m_index(-1)
{
}

KNNGraphQuery::KNNGraphQuery(const KNNGraph* graph, int index) :
    m_graph(graph),
    m_index(index)
{
}

KNNGraphQuery::iterator KNNGraphQuery::begin() const
{
    return m_graph->index_data().begin() + m_index * m_graph->k();
}

KNNGraphQuery::iterator KNNGraphQuery::end() const
{
    return m_graph->index_data().begin() + (m_index+1) * m_graph->k();
}

} // namespace pca
