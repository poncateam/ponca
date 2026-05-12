/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"
#include "../Iterator/neighborGraphKNearestIterator.h"

namespace Ponca
{

#ifndef PARSED_WITH_DOXYGEN
    struct NeighborGraphQueryOutputType : public QueryOutputBase
    {
        using OutputParameter = typename QueryOutputBase::DummyOutputParameter;
    };
#endif

    /*!
     * \brief Extension of the Query class that allows to read the result of a k-nearest neighbors search on the
     * KnnGraph.
     *
     *  Output result of a `KnnGraph::kNearestNeighbors` query request.
     *
     *  \see KnnGraphBase
     */
    template <typename _NeighborGraph>
    class NeighborGraphKNearestQuery
#ifdef PARSED_WITH_DOXYGEN
        : public KNearestIndexQuery<typename NeighborGraph::Traits::IndexType, typename Traits::DataPoint::Scalar>
#else
        // we skip output because we don't need it: k is static, and already stored in the index array
        : public Query<QueryInputIsIndex<typename _NeighborGraph::Traits::IndexType>, NeighborGraphQueryOutputType>
#endif
    {
    public:
        using NeighborGraph = _NeighborGraph;
        using Traits        = typename NeighborGraph::Traits;
        using Iterator = NeighborGraphKNearestIterator<const typename Traits::IndexType*, typename Traits::IndexType>;
#ifdef PARSED_WITH_DOXYGEN
        using QueryType = KNearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>;
#else
        using QueryType = Query<QueryInputIsIndex<typename Traits::IndexType>, NeighborGraphQueryOutputType>;
#endif
        using Self = NeighborGraphKNearestQuery<NeighborGraph>;

    public:
        PONCA_MULTIARCH inline NeighborGraphKNearestQuery(const NeighborGraph* graph, int index)
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
