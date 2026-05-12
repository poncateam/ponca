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
    template <typename _Traits, template <typename> typename BufferType>
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
        using IndexContainerRef = typename Traits::IndexContainerRef; /*!< Ref type to index container */

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

} // namespace Ponca

#undef WRITE_TRAITS
