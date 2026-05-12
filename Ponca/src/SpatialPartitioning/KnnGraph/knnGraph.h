/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./knnGraphTraits.h"

#include "Query/knnGraphKNearestQuery.h"
#include "Query/knnGraphRangeQuery.h"

#include "../KdTree/kdTree.h"
#include "../../Common/Assert.h"

#include <memory>

namespace Ponca
{

    template <typename Traits>
    class KnnGraphBase;

    /*!
     * \brief Public interface for KnnGraph datastructure.
     *
     * Provides default implementation of the KnnGraph
     *
     * \see KnnGraphDefaultTraits for the default trait interface documentation.
     * \see KnnGraphBase for complete API
     */
    template <typename DataPoint>
    using KnnGraph = KnnGraphBase<KnnGraphDefaultTraits<DataPoint>>;

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
    template <typename _Traits, template <typename> typename BufferType>
    class NeighborGraphBase
    {
    public:
#define WRITE_TRAITS                                                                                                    \
    using Traits            = _Traits;                        /*!< Alias to the Traits type                         */  \
    using DataPoint         = typename Traits::DataPoint;     /*!< DataPoint given by user via Traits               */  \
    using Scalar            = typename DataPoint::Scalar;     /*!< Scalar given by user via DataPoint               */  \
    using VectorType        = typename DataPoint::VectorType; /*!< VectorType given by user via DataPoint           */  \
    using IndexType         = typename Traits::IndexType;     /*!< Type used to index points into the PointContainer*/  \
    using PointContainer    = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree   */ \
    using IndexContainer    = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree     */ \
    using IndexContainerRef = typename Traits::IndexContainerRef; /*!< Ref type to index container                */
        WRITE_TRAITS

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

    protected:          // for friends relations
        Buffers m_bufs; ///< Buffers used to store the KnnGraph
    };

    /// \brief Buffer class for StaticKnnGraphBase.
    ///
    /// Extends the NeighborGraphBufferBase with a constant neighborhood size k
    template <typename _Traits>
    struct KnnGraphBuffers : public NeighborGraphBufferBase<_Traits>
    {
        WRITE_TRAITS
        using Base = NeighborGraphBufferBase<_Traits>;

        int k{0};

        PONCA_MULTIARCH inline KnnGraphBuffers() = default;
        PONCA_MULTIARCH inline KnnGraphBuffers(PointContainer _points, const int _k) : Base(_points), k(_k) {}
        PONCA_MULTIARCH inline KnnGraphBuffers(PointContainer _points, typename Traits::IndexContainerRef _indices,
                                               const size_t _points_size, const size_t _indices_size, const int _k)
            : Base(_points, _indices, _points_size, _indices_size), k(_k)
        {
        }
    };

    /*!
     * \brief Customizable base class for KnnGraph datastructure
     *
     * \see Ponca::KnnGraph, Ponca::NeighborGraphBase
     *
     * \tparam _Traits Traits type providing the types and constants used by the KnnGraph. Must have the
     * same interface as the default traits type.
     *
     * \see KnnGraphDefaultTraits for the trait interface documentation.
     *
     */
    template <typename _Traits>
    class StaticKnnGraphBase : public NeighborGraphBase<_Traits, KnnGraphBuffers>
    {
    public:
        WRITE_TRAITS

        using KNearestIndexQuery = KnnGraphKNearestQuery<Traits>;
        using RangeIndexQuery    = KnnGraphRangeQuery<Traits>;
        friend class KnnGraphKNearestQuery<Traits>; /*!< This type must be equal to KnnGraphBase::KNearestIndexQuery
                                                       \see KnnGraphKNearestQuery */
        friend class KnnGraphRangeQuery<Traits>;    /*!< This type must be equal to KnnGraphBase::RangeIndexQuery \see
                                                       KnnGraphRangeQuery */
        using Base    = NeighborGraphBase<Traits, KnnGraphBuffers>;
        using Buffers = typename Base::Buffers;

        PONCA_MULTIARCH inline StaticKnnGraphBase<Traits>(Buffers& _bufs) : Base(_bufs) {}

    protected:
        PONCA_MULTIARCH inline StaticKnnGraphBase(PointContainer _points, const int _k) : Base(Buffers(_points, _k)) {}

        // Query -------------------------------------------------------------------
    public:
        /// \brief Computes a Query object to iterate over the k-nearest neighbors of a point.
        ///
        /// As k was set during the construction of the \ref KnnGraphBase, it doesn't need to be provided.
        ///
        /// The returned object can be reset and reused with the () operator, to compute a new result
        /// (also takes an index as parameter).
        ///
        /// \param index Index of the point that the query evaluates
        /// \return The \ref KNearestIndexQuery mutable object to iterate over the search results.
        PONCA_MULTIARCH inline KNearestIndexQuery kNearestNeighbors(int index) const
        {
            return KNearestIndexQuery(this, index);
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
            return RangeIndexQuery(this, r, index);
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
            return KNearestIndexQuery(this, 0);
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
        PONCA_MULTIARCH inline RangeIndexQuery rangeNeighborsIndexQuery() const { return RangeIndexQuery(this, 0, 0); }

        // Accessors ---------------------------------------------------------------
    public:
        /// \brief Number of neighbor per vertex for a given element (in the KnnGraph, all points have the same number
        /// of neighbors.
        PONCA_MULTIARCH [[nodiscard]] inline int k(int /*index*/ = 0) const { return Base::buffers().k; }
    };

    template <typename _Traits>
    class KnnGraphBase : public StaticKnnGraphBase<_Traits>
    {
    public:
        WRITE_TRAITS
    private:
        using Base    = StaticKnnGraphBase<Traits>;
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
        PONCA_MULTIARCH_HOST inline KnnGraphBase(const KdTreeBase<KdTreeTraits>& _kdtree, const int _k = 6)
            : Base(_kdtree.points(), std::min(_k, _kdtree.sampleCount() - 1))
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

            Base::m_bufs.indices_size = cloudSize * Base::m_bufs.k;
            Base::m_bufs.indices.resize(Base::m_bufs.indices_size, -1);

#pragma omp parallel for shared(_kdtree, cloudSize) default(none)
            for (int i = 0; i < cloudSize; ++i)
            {
                int j = 0;
                for (auto n : _kdtree.kNearestNeighbors(typename KdTreeTraits::IndexType(i),
                                                        typename KdTreeTraits::IndexType(Base::m_bufs.k)))
                {
                    Base::m_bufs.indices[i * Base::m_bufs.k + j] = n;
                    ++j;
                }
            }
        }
    };

} // namespace Ponca

#undef WRITE_TRAITS
