/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "neighborGraphTraits.h"

#include <cstddef> // size_t

#define WRITE_NEIGHBOR_GRAPH_ALIASES                                                                                    \
    using Traits            = _Traits;                    /*!< Alias to the Traits type                          */     \
    using DataPoint         = typename Traits::DataPoint; /*!< DataPoint given by user via Traits                */     \
    using Scalar            = typename DataPoint::Scalar; /*!< Scalar given by user via DataPoint                */     \
    using VectorType        = typename DataPoint::VectorType; /*!< VectorType given by user via DataPoint            */ \
    using IndexType         = typename Traits::IndexType; /*!< Type used to index points into the PointContainer */     \
    using PointContainer    = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree   */ \
    using IndexContainer    = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree     */ \
    using IndexContainerRef = typename Traits::IndexContainerRef; /*!< Ref type to index container */

namespace Ponca
{
    /// \brief Internal structure storing the buffers used by a neighbor graph
    ///
    /// \tparam _Traits Traits of the neighbor graph
    /// \see NeighborGraphDefaultTraits, NeighborGraphPointerTraits
    template <typename _Traits>
    struct NeighborGraphBufferBase
    {
        using Traits         = _Traits;
        using PointContainer = typename Traits::PointContainer;
        using IndexContainer = typename Traits::IndexContainer;

        /// \brief Buffer storing the input points (read only)
        PointContainer points;
        /// \brief Buffer storing the indices associating the input points to the nodes
        IndexContainer indices;

        /// \brief Number of points in the graph
        size_t points_size{0};
        /// \brief Number of connections in the graph
        size_t indices_size{0};

        /// \brief Default constructor, might be deleted depending on the PointContainer type (e.g. const reference)
        PONCA_MULTIARCH inline NeighborGraphBufferBase() = default;

        /// \brief Default constructor setting the point container
        PONCA_MULTIARCH inline NeighborGraphBufferBase(PointContainer _points) : points(_points) {}

        /// \brief Constructor allowing to set all the attributes
        PONCA_MULTIARCH inline NeighborGraphBufferBase(PointContainer _points,
                                                       typename Traits::IndexContainerRef _indices,
                                                       const size_t _points_size, const size_t _indices_size)
            : points(_points), indices(_indices), points_size(_points_size), indices_size(_indices_size)
        {
        }
    };

    /*! \brief Base class for neighbor graphs.
     *
     * This class cannot be used directly: it provides base functionalities (accessors, query generations, buffer
     * storage) and should be inherited to provide:
     *  - procedures to compute the graph,
     *  - methods to access the neighbors of each point, see for instance the class StaticNeighborGraphBase with:
     *  \snippet neighborGraph.h Neighbor Graph Accessors
     *
     *
     * \tparam _Traits Traits type providing the types and constants used by the neighbor graph. Must have the
     * same interface as the default traits types (NeighborGraphDefaultTraits or NeighborGraphPointerTraits).
     * \tparam BufferType Type of buffer used in the Graph. Must inherit NeighborGraphBufferBase and be templated by
     * _Traits
     * \tparam _OneConnectedIndexQuery Type of query used to access direct neighbors
     * \tparam _RangeIndexQuery Type of query used to compute Euclidean range queries restricted on the graph
     */
    template <typename _Traits, template <typename> typename BufferType, typename _OneConnectedIndexQuery,
              typename _RangeIndexQuery>
    class AbstractNeighborGraph
    {
    public:
        WRITE_NEIGHBOR_GRAPH_ALIASES
        using OneConnectedIndexQuery = _OneConnectedIndexQuery;
        using RangeIndexQuery        = _RangeIndexQuery;

        using Buffers = BufferType<Traits>;
        static_assert(std::is_base_of_v<NeighborGraphBufferBase<Traits>, Buffers>,
                      "BufferType must inherit NeighborGraphBufferBase");

    protected:
        /// Accessor used by friend classes (queries) to get const access to the indices whatever the buffer type
        PONCA_MULTIARCH inline const IndexType* getIndexPtr() const { return Traits::getIndexRawPtr(m_bufs.indices); }
        /// Accessor used by friend classes (queries) to get access to the indices whatever the buffer type
        PONCA_MULTIARCH inline IndexType* getIndexPtr() { return Traits::getIndexRawPtr(m_bufs.indices); }

    public:
        // Data --------------------------------------------------------------------

        /*! \brief Constructor that allows the use of prebuilt graph containers.
         *
         * Each internal values of a neighbor graph can be extracted using \ref buffers()
         *
         * \note This constructor can be used to avoid the convertion and building process,
         * which is useful to transfer directly the neighbor graph to the device in CUDA/Sycl.
         *
         * \param _bufs Internal buffers of the AbstractNeighborGraph
         */
        PONCA_MULTIARCH inline AbstractNeighborGraph(const Buffers& _bufs) : m_bufs(_bufs) {}

        //! \brief Get the number of connection edges in the graph
        PONCA_MULTIARCH [[nodiscard]] inline IndexType edgeCount() const { return (IndexType)m_bufs.indices_size; }
        //! \brief Get the number of points
        PONCA_MULTIARCH [[nodiscard]] inline IndexType pointCount() const { return (IndexType)m_bufs.points_size; }
        //! \brief Get the internal point container
        PONCA_MULTIARCH [[nodiscard]] inline PointContainer points() const { return m_bufs.points; };
        //! \brief Get the internal index container
        PONCA_MULTIARCH [[nodiscard]] inline IndexContainer edges() const { return m_bufs.indices; };
        //! \brief Get access to the internal buffer, for instance to prepare GPU binding
        PONCA_MULTIARCH [[nodiscard]] inline const Buffers& buffers() const { return m_bufs; }

        // Query -------------------------------------------------------------------

        /// \brief Provides a Query object to iterate over the vertices that are directly connected to the query point
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index as parameter).
        ///
        /// \param index Index of the point that the query evaluates
        /// \return The \ref OneConnectedIndexQuery mutable object to iterate over the search results.
        PONCA_MULTIARCH [[nodiscard]] inline OneConnectedIndexQuery oneConnectedNeighbors(int index = 0) const
        {
            return OneConnectedIndexQuery(static_cast<const typename OneConnectedIndexQuery::NeighborGraph*>(this),
                                          index);
        }

        /// \brief Provides a Query object to iterate over the neighbors that are inside a given radius.
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index and a radius as parameters).
        ///
        /// \param index Index of the point that the query evaluates
        /// \param r Radius around where to search the neighbors
        /// \return The \ref RangeIndexQuery mutable object to iterate over the search results.
        PONCA_MULTIARCH [[nodiscard]] inline RangeIndexQuery rangeNeighbors(int index, Scalar r) const
        {
            return RangeIndexQuery(static_cast<const typename RangeIndexQuery::NeighborGraph*>(this), r, index);
        }

        /// \brief Convenience function to return an empty mutable range query
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index and a radius as parameters).
        PONCA_MULTIARCH [[nodiscard]] inline RangeIndexQuery rangeNeighborsIndexQuery() const
        {
            return RangeIndexQuery(static_cast<const typename RangeIndexQuery::NeighborGraph*>(this), 0, 0);
        }

    protected:          // for friends relations
        Buffers m_bufs; ///< Buffers used to store the KnnGraph
    };
} // namespace Ponca

