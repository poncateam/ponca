/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "neighborGraphTraits.h"

namespace Ponca
{
    /// \brief Internal structure storing the buffers used by a neighbor graph
    template <typename _Traits>
    struct NeighborGraphBufferBase
    {
        using Traits         = _Traits;
        using PointContainer = typename Traits::PointContainer;
        using IndexContainer = typename Traits::IndexContainer;

        PointContainer points;  ///< Buffer storing the input points (read only)
        IndexContainer indices; ///< Buffer storing the indices associating the input points to the nodes

        size_t points_size{0};
        size_t indices_size{0};

        PONCA_MULTIARCH inline NeighborGraphBufferBase() = default;
        PONCA_MULTIARCH inline NeighborGraphBufferBase(PointContainer _points) : points(_points) {}

        PONCA_MULTIARCH inline NeighborGraphBufferBase(PointContainer _points,
                                                       typename Traits::IndexContainerRef _indices,
                                                       const size_t _points_size, const size_t _indices_size)
            : points(_points), indices(_indices), points_size(_points_size), indices_size(_indices_size)
        {
        }
    };

    /*! \brief Base class for neighbor classes
     *
     * \tparam _Traits Traits type providing the types and constants used by the KnnGraph. Must have the
     * same interface as the default traits type.
     * \tparam BufferType Type of buffer used in the Graph. Must inherit NeighborGraphBufferBase and be templated by
     * _Traits
     */
    template <typename _Traits, template <typename> typename BufferType, typename _KNearestIndexQuery>
    class NeighborGraphBase
    {
    public:
        using Traits         = _Traits;                         /*!< Alias to the Traits type                         */
        using DataPoint      = typename Traits::DataPoint;      /*!< DataPoint given by user via Traits               */
        using Scalar         = typename DataPoint::Scalar;      /*!< Scalar given by user via DataPoint               */
        using VectorType     = typename DataPoint::VectorType;  /*!< VectorType given by user via DataPoint           */
        using IndexType      = typename Traits::IndexType;      /*!< Type used to index points into the PointContainer*/
        using PointContainer = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree  */
        using IndexContainer = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree    */
        using IndexContainerRef  = typename Traits::IndexContainerRef; /*!< Ref type to index container */
        using KNearestIndexQuery = _KNearestIndexQuery;

        using Buffers = BufferType<Traits>;
        static_assert(std::is_base_of_v<NeighborGraphBufferBase<Traits>, Buffers>,
                      "BufferType must inherit NeighborGraphBufferBase");

    protected:
        PONCA_MULTIARCH inline const IndexType* getIndexPtr() const { return Traits::getIndexRawPtr(m_bufs.indices); }
        PONCA_MULTIARCH inline IndexType* getIndexPtr() { return Traits::getIndexRawPtr(m_bufs.indices); }

    public:
        // Data --------------------------------------------------------------------

        /*! \brief Constructor that allows the use of prebuilt KnnGraph containers.
         *
         * Each internal values of a KnnGraph can be extracted using \ref `KnnGraph::buffers()`
         *
         * \note This constructor can be used to avoid the convertion and building process,
         * which is useful to transfer directly the KnnGraph to the device in CUDA.
         *
         * \param _bufs Internal buffers of the KnnGraph
         */
        PONCA_MULTIARCH inline NeighborGraphBase(const Buffers& _bufs) : m_bufs(_bufs) {}

        //! \brief Get the number of indices
        PONCA_MULTIARCH [[nodiscard]] inline IndexType sampleCount() const { return (IndexType)m_bufs.indices_size; }
        //! \brief Get the number of points
        PONCA_MULTIARCH [[nodiscard]] inline IndexType pointCount() const { return (IndexType)m_bufs.points_size; }
        //! \brief Get the internal point container
        PONCA_MULTIARCH [[nodiscard]] inline PointContainer points() const { return m_bufs.points; };
        //! \brief Get the internal index container
        PONCA_MULTIARCH [[nodiscard]] inline IndexContainer samples() const { return m_bufs.indices; };
        //! \brief Get access to the internal buffer, for instance to prepare GPU binding
        PONCA_MULTIARCH [[nodiscard]] inline const Buffers& buffers() const { return m_bufs; }

        // Query -------------------------------------------------------------------

        /// \brief Computes a Query object to iterate over the k-nearest neighbors of a point.
        ///
        /// Here, k is defined by the number of neighbors of the current point in the graph, hence it doesn't need to be
        /// provided.
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index as parameter).
        ///
        /// \param index Index of the point that the query evaluates
        /// \return The \ref KNearestIndexQuery mutable object to iterate over the search results.
        PONCA_MULTIARCH inline KNearestIndexQuery kNearestNeighbors(int index) const
        {
            return KNearestIndexQuery(static_cast<const typename KNearestIndexQuery::NeighborGraph*>(this), index);
        }

        /// \brief Convenience function that provides a k-nearest neighbors Query object.
        ///
        /// Same as `KnnGraphBase::kNearestNeighbors (0)`.
        ///
        /// \warning Unlike `KdTreeBase::kNearestNeighborsIndexQuery`, this function doesn't really return an empty
        /// query. This is due to the fact that the `KnnGraphBase::kNearestNeighbors` query can't set the `k` value to
        /// zero, as it is a value that is managed by the KnnGraphBase structure. Therefore, this function returns the
        /// k-nearest neighbors query made with the evaluation point set to 0.
        ///
        /// \return The \ref KNearestIndexQuery mutable object that can be called with the operator ()
        /// with an index as argument, to fetch the k-nearest neighbors of a point.
        /// \see #kNearestNeighbors
        PONCA_MULTIARCH inline KNearestIndexQuery kNearestNeighborsIndexQuery() const
        {
            return KNearestIndexQuery(static_cast<const typename KNearestIndexQuery::NeighborGraph*>(this), 0);
        }

    protected:          // for friends relations
        Buffers m_bufs; ///< Buffers used to store the KnnGraph
    };

} // namespace Ponca

#undef WRITE_TRAITS
