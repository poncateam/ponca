/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./kdTreeTraits.h"

#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Geometry> // aabb

#include "../../Common/Assert.h"

#include "Query/kdTreeNearestQueries.h"
#include "Query/kdTreeKNearestQueries.h"
#include "Query/kdTreeRangeQueries.h"

namespace Ponca {
template <typename Traits> class KdTreeBase;
template <typename Traits> class KdTreeDenseBase;
template <typename Traits> class KdTreeSparseBase;

/*!
 * \brief Base type for default KdTree implementations
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeBase for complete API
 */
template <typename DataPoint>
using KdTree = KdTreeBase<KdTreeDefaultTraits<DataPoint>>;

/*!
 * \brief Public interface for dense KdTree datastructure.
 *
 * Provides default implementation of the dense KdTree
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeDenseBase for complete API
 */
template <typename DataPoint>
using KdTreeDense = KdTreeDenseBase<KdTreeDefaultTraits<DataPoint>>;

/*!
 * \brief Public interface for sparse KdTree datastructure.
 *
 * Provides default implementation of the sparse KdTree
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeSparseBase for complete API
 */
template <typename DataPoint>
using KdTreeSparse = KdTreeSparseBase<KdTreeDefaultTraits<DataPoint>>;

/*!
 * \brief Customizable base class for KdTree datastructure implementations
 *
 * \see Ponca::KdTreeDense
 * \see Ponca::KdTreeSparse
 *
 * \tparam Traits Traits type providing the types and constants used by the kd-tree. Must have the
 * same interface as the default traits type.
 *
 * \see KdTreeDefaultTraits for the trait interface documentation.
 *
 * \todo Better handle sampling: do not store non-selected points (requires to store original indices)
 */
template <typename Traits>
class KdTreeBase
{
public:
    using DataPoint      = typename Traits::DataPoint; ///< DataPoint given by user via Traits
    using IndexType      = typename Traits::IndexType; ///< Type used to index points into the PointContainer
    using LeafSizeType   = typename Traits::LeafSizeType; ///< Type used to store the size of leaf nodes
    using PointContainer = typename Traits::PointContainer; ///< Container for DataPoint used inside the KdTree
    using IndexContainer = typename Traits::IndexContainer; ///< Container for indices used inside the KdTree
    using NodeIndexType  = typename Traits::NodeIndexType; ///< Type used to index nodes into the NodeContainer
    using NodeType       = typename Traits::NodeType; ///< Type of nodes used inside the KdTree
    using NodeContainer  = typename Traits::NodeContainer; ///< Container for nodes used inside the KdTree

    using Scalar     = typename DataPoint::Scalar; ///< Scalar given by user via DataPoint
    using VectorType = typename DataPoint::VectorType; ///< VectorType given by user via DataPoint
    using AabbType   = typename NodeType::AabbType; ///< Bounding box type given by user via NodeType

    /// \brief The maximum number of nodes that the kd-tree can have.
    static constexpr std::size_t MAX_NODE_COUNT = NodeType::MAX_COUNT;
    /// \brief The maximum number of points that can be stored in the kd-tree.
    static constexpr std::size_t MAX_POINT_COUNT = std::size_t(2) << sizeof(IndexType)*8;

    /// \brief The maximum depth of the kd-tree.
    static constexpr int MAX_DEPTH = Traits::MAX_DEPTH;

    static constexpr bool SUPPORTS_SUBSAMPLING = false;

    static_assert(std::is_same<typename PointContainer::value_type, DataPoint>::value,
        "PointContainer must contain DataPoints");
    
    // Queries use a value of -1 for invalid indices
    static_assert(std::is_signed<IndexType>::value, "Index type must be signed");

    static_assert(std::is_same<typename IndexContainer::value_type, IndexType>::value, "Index type mismatch");
    static_assert(std::is_same<typename NodeContainer::value_type, NodeType>::value, "Node type mismatch");

    static_assert(MAX_DEPTH > 0, "Max depth must be strictly positive");

    // Construction ------------------------------------------------------------
public:
    /// Generate a tree from a custom contained type converted using the specified converter
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling container, transformed to IndexContainer
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename Converter>
    inline void build(PointUserContainer&& points, Converter c);

    /// Convert a custom point container to the KdTree \ref PointContainer using \ref DataPoint default constructor
    struct DefaultConverter
    {
        template <typename Input>
        inline void operator()(Input&& i, PointContainer& o)
        {
            using InputContainer = typename std::remove_reference<Input>::type;
            if constexpr (std::is_same<InputContainer, PointContainer>::value)
                o = std::forward<Input>(i); // Either move or copy
            else
                std::transform(i.cbegin(), i.cend(), std::back_inserter(o),
                               [](const typename InputContainer::value_type &p) -> DataPoint { return DataPoint(p); });
        }
    };

    /// Generate a tree from a custom contained type converted using DefaultConverter
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \param points Input points
    template<typename PointUserContainer>
    inline void build(PointUserContainer&& points)
    {
        build(std::forward<PointUserContainer>(points), DefaultConverter());
    }

    /// Clear tree data
    inline void clear();

    // Accessors ---------------------------------------------------------------
public:
    inline NodeIndexType node_count() const
    {
        return m_nodes.size();
    }

    inline IndexType sample_count() const
    {
        return (IndexType)m_indices.size();
    }

    inline IndexType point_count() const
    {
        return (IndexType)m_points.size();
    }

    inline NodeIndexType leaf_count() const
    {
        return m_leaf_count;
    }

    inline PointContainer& points()
    {
        return m_points;
    };

    inline const PointContainer& points() const
    {
        return m_points;
    };

    inline const NodeContainer& nodes() const
    {
        return m_nodes;
    }

    inline const IndexContainer& samples() const
    {
        return m_indices;
    }

    // Parameters --------------------------------------------------------------
public:
    /// Read leaf min size
    inline LeafSizeType min_cell_size() const
    {
        return m_min_cell_size;
    }

    /// Write leaf min size
    inline void set_min_cell_size(LeafSizeType min_cell_size)
    {
        PONCA_DEBUG_ASSERT(min_cell_size > 0);
        m_min_cell_size = min_cell_size;
    }

    // Index mapping -----------------------------------------------------------
public:
    /// Return the point index associated with the specified sample index
    inline IndexType pointFromSample(IndexType sample_index) const
    {
        return m_indices[sample_index];
    }

    /// Return the \ref DataPoint associated with the specified sample index
    /// \note Convenience function, equivalent to
    /// `point_data()[pointFromSample(sample_index)]`
    inline DataPoint& pointDataFromSample(IndexType sample_index)
    {
        return m_points[pointFromSample(sample_index)];
    }
    
    /// Return the \ref DataPoint associated with the specified sample index
    /// \note Convenience function, equivalent to
    /// `point_data()[pointFromSample(sample_index)]`
    inline const DataPoint& pointDataFromSample(IndexType sample_index) const
    {
        return m_points[pointFromSample(sample_index)];
    }

    // Query -------------------------------------------------------------------
public :
    KdTreeKNearestPointQuery<Traits> k_nearest_neighbors(const VectorType& point, IndexType k) const
    {
        return KdTreeKNearestPointQuery<Traits>(this, k, point);
    }

    KdTreeKNearestIndexQuery<Traits> k_nearest_neighbors(IndexType index, IndexType k) const
    {
        return KdTreeKNearestIndexQuery<Traits>(this, k, index);
    }

    KdTreeNearestPointQuery<Traits> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<Traits>(this, point);
    }

    KdTreeNearestIndexQuery<Traits> nearest_neighbor(IndexType index) const
    {
        return KdTreeNearestIndexQuery<Traits>(this, index);
    }

    KdTreeRangePointQuery<Traits> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<Traits>(this, r, point);
    }

    KdTreeRangeIndexQuery<Traits> range_neighbors(IndexType index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<Traits>(this, r, index);
    }
    
    // Utilities ---------------------------------------------------------------
public:
    inline bool valid() const;
    inline void print(std::ostream& os, bool verbose = false) const;

    // Data --------------------------------------------------------------------
protected:
    PointContainer m_points;
    NodeContainer m_nodes;
    IndexContainer m_indices;

    LeafSizeType m_min_cell_size {64}; ///< Minimal number of points per leaf
    NodeIndexType m_leaf_count {0}; ///< Number of leaves in the Kdtree (computed during construction)

    // Internal ----------------------------------------------------------------
protected:
    inline KdTreeBase() = default;

    /// Generate a tree sampled from a custom contained type converted using a `Converter`
    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam Converter
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
    inline void buildWithSampling(PointUserContainer&& points,
                                  IndexUserContainer sampling,
                                  Converter c);

    /// Generate a tree sampled from a custom contained type converted using a \ref KdTreeBase::DefaultConverter
    /// \tparam PointUserContainer Input points, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Samples used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline void buildWithSampling(PointUserContainer&& points,
                                  IndexUserContainer sampling)
    {
        buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling), DefaultConverter());
    }

private:
    inline void build_rec(NodeIndexType node_id, IndexType start, IndexType end, int level);
    inline IndexType partition(IndexType start, IndexType end, int dim, Scalar value);
};

/*!
 * \brief Customizable base class for dense KdTree datastructure
 *
 * This version of the KdTree does not support subsampling. For an
 * implementation that supports subsampling, see \ref KdTreeSparseBase.
 *
 * \see Ponca::KdTreeDense
 * \see Ponca::KdTreeSparse
 *
 * \tparam Traits Traits type providing the types and constants used by the kd-tree. Must have the
 * same interface as the default traits type.
 *
 * \see KdTreeDefaultTraits for the trait interface documentation.
 */
template <typename Traits>
class KdTreeDenseBase : public KdTreeBase<Traits>
{
private:
    using Base = KdTreeBase<Traits>;

public:
    /// Default constructor creating an empty tree
    /// \see build
    KdTreeDenseBase() = default;

    /// Constructor generating a tree from a custom contained type converted using a \ref KdTreeBase::DefaultConverter
    template<typename PointUserContainer>
    inline explicit KdTreeDenseBase(PointUserContainer&& points)
        : Base()
    {
        this->build(std::forward<PointUserContainer>(points));
    }
};

/*!
 * \brief Customizable base class for KdTreeSparse datastructure
 *
 * This version of the KdTree supports construction using a subset of samples.
 *
 * \see buildWithSampling
 * \see Ponca::KdTreeDense
 * \see Ponca::KdTreeSparse
 *
 * \tparam Traits Traits type providing the types and constants used by the kd-tree. Must have the
 * same interface as the default traits type.
 *
 * \see KdTreeDefaultTraits for the trait interface documentation.
 */
template <typename Traits>
class KdTreeSparseBase : public KdTreeBase<Traits>
{
private:
    using Base = KdTreeBase<Traits>;

public:
    static constexpr bool SUPPORTS_SUBSAMPLING = false;

    /// Default constructor creating an empty tree
    /// \see build
    KdTreeSparseBase() = default;

    /// Constructor generating a tree from a custom contained type converted using a \ref KdTreeBase::DefaultConverter
    template<typename PointUserContainer>
    inline explicit KdTreeSparseBase(PointUserContainer&& points)
        : Base()
    {
        this->build(std::forward<PointUserContainer>(points));
    }

    /// Constructor generating a tree sampled from a custom contained type converted using a \ref KdTreeBase::DefaultConverter
    /// \tparam PointUserContainer Input points, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param point Input points
    /// \param sampling Samples used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTreeSparseBase(PointUserContainer&& points, IndexUserContainer sampling)
        : Base()
    {
        this->buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling));
    }

    using Base::buildWithSampling;
};

#include "./kdTree.hpp"
} // namespace Ponca

template <typename Traits>
std::ostream& operator<<(std::ostream& os, const Ponca::KdTreeBase<Traits>& kdtree)
{
    kdtree.print(os);
    return os;
}
