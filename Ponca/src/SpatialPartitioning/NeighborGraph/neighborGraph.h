/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "neighborGraphTraits.h"
#include "Query/neighborGraphKNearestQuery.h"
#include "Query/neighborGraphRangeQuery.h"

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
    template <typename _Traits, template <typename> typename BufferType, typename _KNearestIndexQuery,
              typename _RangeIndexQuery>
    class AbstractNeighborGraph
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
        using RangeIndexQuery    = _RangeIndexQuery;

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
        PONCA_MULTIARCH inline AbstractNeighborGraph(const Buffers& _bufs) : m_bufs(_bufs) {}

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

        /// \brief Computes a Query object to iterate over the neighbors that are inside a given radius.
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index and a radius as parameters).
        ///
        /// \param index Index of the point that the query evaluates
        /// \param r Radius around where to search the neighbors
        /// \return The \ref RangeIndexQuery mutable object to iterate over the search results.
        PONCA_MULTIARCH inline RangeIndexQuery rangeNeighbors(int index, Scalar r) const
        {
            return RangeIndexQuery(static_cast<const typename RangeIndexQuery::NeighborGraph*>(this), r, index);
        }

        /// \brief Convenience function that provides an empty range neighbors Query object.
        ///
        /// The returned object can be called with the arguments `(i, r)` to fetch the neighbors
        /// that are in range `r` of the point of index `i`.
        ///
        /// Same as `KnnGraphBase::rangeNeighbors (0, 0)`.
        ///
        /// \return The empty \ref KNearestIndexQuery mutable object to iterate over the search results.
        /// \see #rangeNeighbors
        PONCA_MULTIARCH inline RangeIndexQuery rangeNeighborsIndexQuery() const
        {
            return RangeIndexQuery(static_cast<const typename RangeIndexQuery::NeighborGraph*>(this), 0, 0);
        }

    protected:          // for friends relations
        Buffers m_bufs; ///< Buffers used to store the KnnGraph
    };

    /// \brief Buffer class for StaticNeighborGraphBase.
    ///
    /// Extends the NeighborGraphBufferBase with a buffer storing the start/end indices of each vertex neighborhood.
    /// For the vertex i, the range is accessed as follows:
    ///   - start(i) = ranges[i]
    ///   - end(i) = ranges[i+1]
    /// \see StaticNeighborGraphBase::beginId, StaticNeighborGraphBase::endId
    template <typename _Traits>
    struct NeighborGraphBuffer : public NeighborGraphBufferBase<_Traits>
    {
#define WRITE_TRAITS                                                                                                    \
    using Traits            = _Traits;                    /*!< Alias to the Traits type                          */     \
    using DataPoint         = typename Traits::DataPoint; /*!< DataPoint given by user via Traits                */     \
    using Scalar            = typename DataPoint::Scalar; /*!< Scalar given by user via DataPoint                */     \
    using VectorType        = typename DataPoint::VectorType; /*!< VectorType given by user via DataPoint            */ \
    using IndexType         = typename Traits::IndexType; /*!< Type used to index points into the PointContainer */     \
    using PointContainer    = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree   */ \
    using IndexContainer    = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree     */ \
    using IndexContainerRef = typename Traits::IndexContainerRef; /*!< Ref type to index container */
        WRITE_TRAITS
        using Base = NeighborGraphBufferBase<_Traits>;

        IndexContainer ranges;

        PONCA_MULTIARCH inline NeighborGraphBuffer() = default;
        PONCA_MULTIARCH inline NeighborGraphBuffer(PointContainer _points) : Base(_points) {}
        PONCA_MULTIARCH inline NeighborGraphBuffer(PointContainer _points, typename Traits::IndexContainerRef _indices,
                                                   const size_t _points_size, const size_t _indices_size,
                                                   typename Traits::IndexContainerRef _ranges)
            : Base(_points, _indices, _points_size, _indices_size), ranges(_ranges)
        {
        }
    };

    /*!
     * \brief Static generic neighbor graph (does not include build functions)
     *
     * In this generic version, the graph holds an index container that stores the begin and end indices of each
     * point's neighborhood
     *
     * \see Ponca::StaticKnnGraphBase
     *
     * \tparam _Traits Traits type providing the types and constants used by the KnnGraph. Must have the
     * same interface as the default traits type.
     *
     * \see NeighborGraphDefaultTraits for the trait interface documentation.
     *
     */
    template <typename _Traits>
    class StaticNeighborGraphBase
        : public AbstractNeighborGraph<_Traits, NeighborGraphBuffer,
                                       NeighborGraphKNearestQuery<StaticNeighborGraphBase<_Traits>>,
                                       NeighborGraphRangeQuery<StaticNeighborGraphBase<_Traits>>>
    {
    public:
        using Traits         = _Traits;                         /*!< Alias to the Traits type                         */
        using DataPoint      = typename Traits::DataPoint;      /*!< DataPoint given by user via Traits               */
        using Scalar         = typename DataPoint::Scalar;      /*!< Scalar given by user via DataPoint               */
        using VectorType     = typename DataPoint::VectorType;  /*!< VectorType given by user via DataPoint           */
        using IndexType      = typename Traits::IndexType;      /*!< Type used to index points into the PointContainer*/
        using PointContainer = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree  */
        using IndexContainer = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree    */
        using IndexContainerRef = typename Traits::IndexContainerRef; /*!< Ref type to index container */

        friend class NeighborGraphKNearestQuery<StaticNeighborGraphBase<Traits>>; /*!< This type must be equal to
                                                       KnnGraphBase::KNearestIndexQuery
                                                       \see NeighborGraphKNearestQuery */
        friend class NeighborGraphRangeQuery<StaticNeighborGraphBase<Traits>>;    /*!< This type must be equal to
                                                       KnnGraphBase::RangeIndexQuery \see   NeighborGraphRangeQuery */
        using Base    = AbstractNeighborGraph<Traits, NeighborGraphBuffer,
                                              NeighborGraphKNearestQuery<StaticNeighborGraphBase<Traits>>,
                                              NeighborGraphRangeQuery<StaticNeighborGraphBase<Traits>>>;
        using Buffers = typename Base::Buffers;

        PONCA_MULTIARCH inline StaticNeighborGraphBase(const Buffers& _bufs) : Base(_bufs) {}

        // Accessors ---------------------------------------------------------------
    public:
        /// \brief Number of neighbor per vertex for a given element (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int k(int vId = 0) const { return endId(vId) * beginId(vId); }
        /// Index of the beginning of the neighborhood range (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int beginId(int vId) const { return Base::buffers().ranges[vId]; }
        /// Index of the end of the neighborhood range (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int endId(int vId) const { return Base::buffers().ranges[vId + 1]; }
    };

    template <typename _Traits>
    class NeighborGraphBase : public StaticNeighborGraphBase<_Traits>
    {
    public:
        WRITE_TRAITS
    private:
        using Base    = StaticNeighborGraphBase<Traits>;
        using Buffers = typename Base::Buffers;
        // knnGraph ----------------------------------------------------------------
    public:
        /// \brief Build a KnnGraph from a KdTreeDense
        ///
        /// \param _kdtree Reference to the KdTree
        /// \param _k Number of requested neighbors. Might be reduced if k is larger than the kdtree size - 1
        ///          (query point is not included in query output, thus -1)
        ///
        /// \warning Stores a const reference to kdtree.point_data()
        /// \warning KdTreeTraits compatibility is checked with static assertion
        template <typename KdTreeTraits>
        PONCA_MULTIARCH_HOST inline NeighborGraphBase(const KdTreeBase<KdTreeTraits>& _kdtree, const Scalar range)
            : Base(_kdtree.points())
        {
            Base::m_bufs.points_size = _kdtree.pointCount();
#define CHECK_TRAITS_TYPENAME_COMPAT(A, B)                                                            \
    static_assert(std::is_same_v<A, B> || std::is_convertible_v<A, B> || std::is_convertible_v<B, A>, \
                  "KdTreeTraits::DataPoint is not equal to Traits::DataPoint");

            static_assert(std::is_same_v<typename Traits::DataPoint, typename KdTreeTraits::DataPoint>,
                          "KdTreeTraits::DataPoint is not equal to Traits::DataPoint");

            CHECK_TRAITS_TYPENAME_COMPAT(typename Traits::PointContainer, typename KdTreeTraits::PointContainer)
            CHECK_TRAITS_TYPENAME_COMPAT(typename Traits::IndexContainer, typename KdTreeTraits::IndexContainer)

#undef CHECK_TRAITS_TYPENAME_COMPAT

            // We need to account for the entire point set, irrespectively of the sampling. This is because the kdtree
            // (kNearestNeighbors) return ids of the entire point set, not it sub-sampled list of ids.
            // \fixme Update API to properly handle kdtree subsampling
            const int cloudSize = _kdtree.pointCount();
            {
                const int samplesSize = _kdtree.sampleCount();
                PONCA_ASSERT(cloudSize == samplesSize);
            }

            Base::m_bufs.ranges.resize(cloudSize + 1); // we need one more index (see StaticNeighborGraphBase::endId)
            Base::m_bufs.ranges[0] = 0;                // first element starts at 0;

            for (int i = 0; i < cloudSize; ++i)
            {
                for (auto n : _kdtree.rangeNeighbors(typename KdTreeTraits::IndexType(i), range))
                {
                    Base::m_bufs.indices.push_back(n);
                }
                Base::m_bufs.ranges[i + 1] = Base::m_bufs.indices.size();
            }
            Base::m_bufs.indices_size = Base::m_bufs.indices.size();
        }
    };

    /*!
     * \brief Public interface for KnnGraph datastructure.
     *
     * Provides default implementation of the KnnGraph
     *
     * \see NeighborGraphDefaultTraits for the default trait interface documentation.
     * \see KnnGraphBase for complete API
     */
    template <typename DataPoint>
    using NeighborGraph = NeighborGraphBase<NeighborGraphDefaultTraits<DataPoint>>;

} // namespace Ponca

#undef WRITE_TRAITS
