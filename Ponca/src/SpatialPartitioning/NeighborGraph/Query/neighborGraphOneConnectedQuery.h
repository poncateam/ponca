/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"
#include "../Iterator/neighborGraphOneConnectedIterator.h"

namespace Ponca
{
    /// \brief Output type of the NeighborGraphOneConnectedQuery
    struct NeighborGraphQueryOutputType : public QueryOutputBase
    {
        using OutputParameter = typename QueryOutputBase::DummyOutputParameter;
    };

    /*!
     * \brief Extension of the Query class that allows to read the neighbors that are directly connected to the query
     * point in the neighbor graph.
     *
     *  Output result of a `KnnGraph::oneConnectedNeighbors` query request.
     *
     *  \see KnnGraphBase
     */
    template <typename _NeighborGraph>
    class NeighborGraphOneConnectedQuery
        : public Query<QueryInputIsIndex<typename _NeighborGraph::Traits::IndexType>, NeighborGraphQueryOutputType>
    {
    public:
        using NeighborGraph = _NeighborGraph;
        using Traits        = typename NeighborGraph::Traits;
        using Iterator =
            NeighborGraphOneConnectedIterator<const typename Traits::IndexType*, typename Traits::IndexType>;
        using QueryType = Query<QueryInputIsIndex<typename Traits::IndexType>, NeighborGraphQueryOutputType>;

        using Self = NeighborGraphOneConnectedQuery<NeighborGraph>;

    public:
        PONCA_MULTIARCH inline NeighborGraphOneConnectedQuery(const NeighborGraph* graph, int index)
            : QueryType(index), m_graph(graph)
        {
        }
        /// \brief Call the k-nearest neighbors query with new input parameter.
        PONCA_MULTIARCH inline Self& operator()(int index) { return QueryType::template operator()<Self>(index); }

        /// \brief Returns an iterator to the beginning of the k-nearest neighbors query.
        PONCA_MULTIARCH [[nodiscard]] inline Iterator begin() const
        {
            return Iterator(m_graph->getIndexPtr()) + m_graph->beginId(QueryType::input());
        }

        /// \brief Returns an iterator to the end of the k-nearest neighbors query.
        PONCA_MULTIARCH [[nodiscard]] inline Iterator end() const
        {
            return Iterator(m_graph->getIndexPtr()) + m_graph->endId(QueryType::input());
        }

    protected:
        const NeighborGraph* m_graph{nullptr};
    };

} // namespace Ponca
