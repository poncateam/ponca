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

    /*!
     * \brief Customizable base class for KnnGraph datastructure
     *
     * \see Ponca::KnGraph
     *
     * \tparam Traits Traits type providing the types and constants used by the KnnGraph. Must have the
     * same interface as the default traits type.
     *
     * \see KnnGraphDefaultTraits for the trait interface documentation.
     *
     */
    template <typename Traits>
    class StaticKnnGraphBase
    {
    public:
#define WRITE_TRAITS                                                                                                 \
    using DataPoint      = typename Traits::DataPoint;      /*!< DataPoint given by user via Traits               */ \
    using Scalar         = typename DataPoint::Scalar;      /*!< Scalar given by user via DataPoint               */ \
    using VectorType     = typename DataPoint::VectorType;  /*!< VectorType given by user via DataPoint           */ \
    using IndexType      = typename Traits::IndexType;      /*!< Type used to index points into the PointContainer*/ \
    using PointContainer = typename Traits::PointContainer; /*!< Container for DataPoint used inside the KdTree   */ \
    using IndexContainer = typename Traits::IndexContainer; /*!< Container for indices used inside the KdTree     */
        WRITE_TRAITS

        using KNearestIndexQuery = KnnGraphKNearestQuery<Traits>;
        using RangeIndexQuery    = KnnGraphRangeQuery<Traits>;
        friend class KnnGraphKNearestQuery<Traits>; /*!< This type must be equal to KnnGraphBase::KNearestIndexQuery*/
        friend class KnnGraphRangeQuery<Traits>;    /*!< This type must be equal to KnnGraphBase::RangeIndexQuery   */

        /// \brief Internal structure storing all the buffers used by the KdTree
        struct Buffers
        {
            PointContainer points;  ///< Buffer storing the input points (read only)
            IndexContainer indices; ///< Buffer storing the indices associating the input points to the nodes

            size_t points_size{0};
            size_t indices_size{0};

            PONCA_MULTIARCH inline Buffers() = default;

            PONCA_MULTIARCH inline Buffers(PointContainer _points, IndexContainer _indices, const size_t _points_size,
                                           const size_t _indices_size)
                : points(_points), indices(_indices), points_size(_points_size), indices_size(_indices_size)
            {
            }
        };

    protected:
        PONCA_MULTIARCH inline StaticKnnGraphBase(int _k) : m_bufs(), m_k(_k){};

    public:
        /*! \brief Constructor that allows the use of prebuilt KnnGraph containers.
         *
         * Each internal values of a KnnGraph can be extracted using \ref `KnnGraph::buffers()`
         *
         * \note This constructor can be used to avoid the convertion and building process,
         * which is useful to transfer directly the KnnGraph to the device in CUDA.
         *
         * \param buf Internal buffers of the KnnGraph
         */
        PONCA_MULTIARCH inline StaticKnnGraphBase(Buffers& _bufs, int _k) : m_bufs(_bufs), m_k(_k) {}

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
        /// \brief Number of neighbor per vertex
        PONCA_MULTIARCH [[nodiscard]] inline int k() const { return m_k; }
        /// \brief Number of vertices in the neighborhood graph
        PONCA_MULTIARCH [[nodiscard]] inline size_t size() const
        {
            return m_bufs.indices_size / static_cast<size_t>(m_k);
        }
        //! \brief Get the number of indices
        PONCA_MULTIARCH [[nodiscard]] inline IndexType sampleCount() const { return (IndexType)m_bufs.indices_size; }
        //! \brief Get the number of points
        PONCA_MULTIARCH [[nodiscard]] inline IndexType pointCount() const { return (IndexType)m_bufs.points_size; }
        //! \brief Get the internal point container
        PONCA_MULTIARCH [[nodiscard]] inline const PointContainer& points() const { return m_bufs.points; };
        //! \brief Get the internal indice container
        PONCA_MULTIARCH [[nodiscard]] inline const IndexContainer& samples() const { return m_bufs.indices; };
        //! \brief Get access to the internal buffer, for instance to prepare GPU binding
        PONCA_MULTIARCH [[nodiscard]] inline const Buffers& buffers() const { return m_bufs; }

        // Data --------------------------------------------------------------------
    protected:          // for friends relations
        Buffers m_bufs; ///< Buffers used to store the KnnGraph
        const int m_k;
    };

    template <typename Traits>
    class KnnGraphBase : public StaticKnnGraphBase<Traits>
    {
    public:
        WRITE_TRAITS
    private:
        using Base = StaticKnnGraphBase<Traits>;
        // knnGraph ----------------------------------------------------------------
    public:
        /// \brief Build a KnnGraph from a KdTreeDense
        ///
        /// \param k Number of requested neighbors. Might be reduced if k is larger than the kdtree size - 1
        ///          (query point is not included in query output, thus -1)
        ///
        /// \warning Stores a const reference to kdtree.point_data()
        /// \warning KdTreeTraits compatibility is checked with static assertion
        template <typename KdTreeTraits>
        PONCA_MULTIARCH inline KnnGraphBase(const KdTreeBase<KdTreeTraits>& _kdtree, const int _k = 6)
            : Base(std::min(_k, _kdtree.sampleCount() - 1))
        {
            Base::m_bufs.points_size = _kdtree.points().size();
            Base::m_bufs.points      = std::move(_kdtree.points());
            static_assert(std::is_same_v<typename Traits::DataPoint, typename KdTreeTraits::DataPoint>,
                          "KdTreeTraits::DataPoint is not equal to Traits::DataPoint");
            static_assert(std::is_same_v<typename Traits::PointContainer, typename KdTreeTraits::PointContainer>,
                          "KdTreeTraits::PointContainer is not equal to Traits::PointContainer");
            static_assert(std::is_same_v<typename Traits::IndexContainer, typename KdTreeTraits::IndexContainer>,
                          "KdTreeTraits::IndexContainer is not equal to Traits::IndexContainer");

            // We need to account for the entire point set, irrespectively of the sampling. This is because the kdtree
            // (kNearestNeighbors) return ids of the entire point set, not it sub-sampled list of ids.
            // \fixme Update API to properly handle kdtree subsampling
            const int cloudSize = _kdtree.pointCount();
            {
                const int samplesSize = _kdtree.sampleCount();
                eigen_assert(cloudSize == samplesSize);
            }

            Base::m_bufs.indices_size = cloudSize * Base::m_k;
            Base::m_bufs.indices.resize(Base::m_bufs.indices_size, -1);

#pragma omp parallel for shared(_kdtree, cloudSize) default(none)
            for (int i = 0; i < cloudSize; ++i)
            {
                int j = 0;
                for (auto n : _kdtree.kNearestNeighbors(typename KdTreeTraits::IndexType(i),
                                                        typename KdTreeTraits::IndexType(Base::m_k)))
                {
                    Base::m_bufs.indices[i * Base::m_k + j] = n;
                    ++j;
                }
            }
        }
    };

} // namespace Ponca

#undef WRITE_TRAITS
