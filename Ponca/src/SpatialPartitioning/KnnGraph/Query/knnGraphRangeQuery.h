/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"
#include "../Iterator/knnGraphRangeIterator.h"
#include "../../../Common/Containers/stack.h"

namespace Ponca
{
    template <typename Traits>
    class StaticKnnGraphBase;

    /*! \brief Extension of the Query class that allows to read the result of a range neighbor search on the KnnGraph.
     *
     * Output result of a `KnnGraph::rangeNeighbors` query request. \see StaticKnnGraphBase
     *
     * Contrary to the `KdTreeRangeQuery` this `RangeIndexQuery` needs to keep track of all the already visited
     * neighbors using a Set-like Data structure.
     *
     * \tparam IndexSet An index set type used to avoid the already visited neighbor in our search algorithm. This
     * structure is required by the `KnnGraphRangeQuery::advance` method.
     *
     * For the `IndexSet`, Ponca provides two possible choices :
     *
     * - (Default) `HashSet<Traits::MAX_RANGE_NEIGHBORS_SIZE>` : Stores the index in a HashMap-like structure.
     * The Best case complexity for insertion and search is O(1) and worst case is O(N), depending on the given dataset
     * and on the chosen hashing function (Sparser hashing results will lead to a reduce look-up time).
     * \see HashSet
     *
     * - `BitSet<MAX_POINT_CLOUD_SIZE>` : Stores the index in a set of bits by allocating a single bit to each possible
     * indices. This structure provides better search and insert complexity at the cost of more memory use (trivial
     * insertion and search, with O(1) complexity). Since we need to allocate a single memory bit for each possible
     * indices, to indicate if the index was stored or not, the memory use of this BitSet type is dependant of the total
     * size of the point cloud.
     * Extensively big point clouds can therefore provoke memory allocation problems if we have a memory limit (e.g. if
     * we are instantiating the `KnnGraphRangeQuery` inside local memory, in a CUDA kernel).
     * \see BitSet for more detailed information about the memory usage.
     */
    template <typename Traits, typename IndexSet>
    class KnnGraphRangeQuery : public RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
    {
    protected:
        using QueryType = RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>;
        friend class KnnGraphRangeIterator<Traits>; // This type must be equal to KnnGraphRangeQuery::Iterator

    public:
        using DataPoint  = typename Traits::DataPoint;
        using IndexType  = typename Traits::IndexType;
        using Scalar     = typename DataPoint::Scalar;
        using VectorType = typename DataPoint::VectorType;
        using Iterator   = KnnGraphRangeIterator<Traits>;
        using Self       = KnnGraphRangeQuery<Traits>;

    public:
        PONCA_MULTIARCH inline KnnGraphRangeQuery(const StaticKnnGraphBase<Traits>* graph, Scalar radius, int index)
            : QueryType(radius, index), m_graph(graph), m_flag(), m_stack()
        {
        }

        /// \brief Call the range neighbors query with new input and radius parameters.
        PONCA_MULTIARCH inline Self& operator()(int index, Scalar radius)
        {
            return QueryType::template operator()<Self>(index, radius);
        }

        /// \brief Call the range neighbors query with new input parameter.
        PONCA_MULTIARCH inline Self& operator()(int index) { return QueryType::template operator()<Self>(index); }

        /// \brief Returns an iterator to the beginning of the range neighbors query.
        PONCA_MULTIARCH inline Iterator begin()
        {
            QueryType::reset();
            Iterator it(this);
            this->initialize(it);
            this->advance(it);
            return it;
        }

        /// \brief Returns an iterator to the end of the range neighbors query.
        PONCA_MULTIARCH inline Iterator end() { return Iterator(this, static_cast<int>(m_graph->size())); }

    protected:
        PONCA_MULTIARCH inline void initialize(Iterator& iterator)
        {
            m_flag.clear();
            m_flag.insert(QueryType::input());

            PONCA_DEBUG_ASSERT(m_stack.empty());
            m_stack.push(QueryType::input());

            iterator.m_index = -1;
        }

        /*! \brief Helper function for the KnnGraphRangeIterator that advances the range neighbors search using the
         * k-nearest neighbors known by the KnnGraph
         *
         * \param iterator The KnnGraphRangeIterator from where the advance request is made
         * \see KnnGraphRangeIterator
         */
        PONCA_MULTIARCH inline void advance(Iterator& iterator)
        {
            const auto& points = m_graph->points();
            const auto& point  = points[QueryType::input()].pos();

            if (!(iterator != end()))
                return;

            if (m_stack.empty())
            {
                iterator = end();
            }
            else
            {
                int idx_current = m_stack.top();
                m_stack.pop();

                PONCA_DEBUG_ASSERT((point - points[idx_current].pos()).squaredNorm() < QueryType::squaredRadius());

                iterator.m_index = idx_current;

                for (int idx_nei : m_graph->kNearestNeighbors(idx_current))
                {
                    PONCA_DEBUG_ASSERT(idx_nei >= 0);
                    Scalar d  = (point - points[idx_nei].pos()).squaredNorm();
                    Scalar th = QueryType::descentDistanceThreshold();

                    if ((point - points[idx_nei].pos()).squaredNorm() < QueryType::descentDistanceThreshold() &&
                        m_flag.insert(idx_nei))
                    {
                        m_stack.push(idx_nei);
                    }
                }
                if (iterator.m_index == QueryType::input())
                    advance(iterator); // query is not included in returned set
            }
        }

    protected:
        const StaticKnnGraphBase<Traits>* m_graph{nullptr};
        IndexSet m_flag;                                      ///< Stores every visited neighbor ids
        Stack<int, Traits::MAX_RANGE_NEIGHBORS_SIZE> m_stack; ///< Holds the next ids the Query should visit
    };

} // namespace Ponca
