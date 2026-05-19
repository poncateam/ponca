/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "abstractNeighborGraph.h"
#include "Query/neighborGraphOneConnectedQuery.h"
#include "Query/neighborGraphRangeQuery.h"

#include "../KdTree/kdTree.h"

namespace Ponca
{
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
        WRITE_NEIGHBOR_GRAPH_ALIASES
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
     * In this generic version, the graph holds an index container that stores the begin/end indices of each
     * point's neighborhood
     *
     * \tparam _Traits Traits type providing the types and constants used by the neighbor graph. Must have the
     * same interface as the default traits type (e.g., NeighborGraphDefaultTraits or NeighborGraphPointerTraits).
     *
     */
    template <typename _Traits>
    class StaticNeighborGraphBase
        : public AbstractNeighborGraph<_Traits, NeighborGraphBuffer,
                                       NeighborGraphOneConnectedQuery<StaticNeighborGraphBase<_Traits>>,
                                       NeighborGraphRangeQuery<StaticNeighborGraphBase<_Traits>>>
    {
    public:
        WRITE_NEIGHBOR_GRAPH_ALIASES

        /// This type must be equal to AbstractNeighborGraph::OneConnectedIndexQuery \see NeighborGraphKNearestQuery
        friend class NeighborGraphOneConnectedQuery<StaticNeighborGraphBase<Traits>>;
        /// This type must be equal to AbstractNeighborGraph::RangeIndexQuery \see   NeighborGraphRangeQuery
        friend class NeighborGraphRangeQuery<StaticNeighborGraphBase<Traits>>;

        using Base    = AbstractNeighborGraph<Traits, NeighborGraphBuffer,
                                              NeighborGraphOneConnectedQuery<StaticNeighborGraphBase<Traits>>,
                                              NeighborGraphRangeQuery<StaticNeighborGraphBase<Traits>>>;
        using Buffers = typename Base::Buffers;

        PONCA_MULTIARCH inline StaticNeighborGraphBase(const Buffers& _bufs) : Base(_bufs) {}

        //! [Neighbor Graph Accessors]
        /// \brief Number of neighbor per vertex for a given element (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int k(int vId = 0) const { return endId(vId) * beginId(vId); }
        /// Index of the beginning of the neighborhood range (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int beginId(int vId) const { return Base::buffers().ranges[vId]; }
        /// Index of the end of the neighborhood range (see storing convention in NeighborGraphBuffer)
        PONCA_MULTIARCH [[nodiscard]] inline int endId(int vId) const { return Base::buffers().ranges[vId + 1]; }
        //! [Neighbor Graph Accessors]
    };

    template <typename _Traits>
    class NeighborGraphBase : public StaticNeighborGraphBase<_Traits>
    {
    public:
        WRITE_NEIGHBOR_GRAPH_ALIASES

    private:
        using Base    = StaticNeighborGraphBase<Traits>;
        using Buffers = typename Base::Buffers;

    public:
        /// \brief Build a neighbor graph by connecting all points using Euclidean range queries: two points
        /// \f[\mathbf{p,q}\f] are connected if \f[\mathbf{p}-\mathbf{q}<\text{range}\f]
        ///
        /// \note Each point might have a different number of neighbors
        /// \note Empty neighborhood are not checked
        ///
        /// \param _kdtree Reference to the KdTree
        /// \param range Distance threshold used to connect points
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
     * \brief Public interface for the NeighborGraph datastructure.
     *
     * Provides default implementation of the NeighborGraph
     *
     * \see NeighborGraphDefaultTraits for the default trait interface documentation.
     * \see NeighborGraphBase, AbstractNeighborGraph for complete API
     */
    template <typename DataPoint>
    using NeighborGraph = NeighborGraphBase<NeighborGraphDefaultTraits<DataPoint>>;

} // namespace Ponca

#undef WRITE_TRAITS
