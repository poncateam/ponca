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
 * \brief Abstract KdTree type with KdTreeDefaultTraits
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeBase for complete API
 *
 * \warning It is not possible to create instances of type KdTree. This type must only be used to store pointers
 * to KdTreeDense or KdTreeSparse objects, e.g. the declaration
 * \snippet examples/cpp/ponca_neighbor_search.cpp KdTree pointer usage
 * \snippet examples/cpp/ponca_neighbor_search.cpp KdTree assign sparse
 * \snippet examples/cpp/ponca_neighbor_search.cpp KdTree assign dense
 */
#ifdef PARSED_WITH_DOXYGEN
/// [KdTree type definition]
template <typename DataPoint>
struct KdTree : public Ponca::KdTreeBase<KdTreeDefaultTraits<DataPoint>>{};
/// [KdTree type definition]
#else
template <typename DataPoint>
using KdTree = KdTreeBase<KdTreeDefaultTraits<DataPoint>>; // prefer alias to avoid redefining methods
#endif


/*!
 * \brief Public interface for dense KdTree datastructure.
 *
 * Provides default implementation of the dense KdTree
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeDenseBase for complete API
 */
#ifdef PARSED_WITH_DOXYGEN
/// [KdTreeDense type definition]
template <typename DataPoint>
struct KdTreeDense : public Ponca::KdTreeDenseBase<KdTreeDefaultTraits<DataPoint>>{};
/// [KdTreeDense type definition]
#else
template <typename DataPoint>
using KdTreeDense = KdTreeDenseBase<KdTreeDefaultTraits<DataPoint>>; // prefer alias to avoid redefining methods
#endif

/*!
 * \brief Public interface for sparse KdTree datastructure.
 *
 * Provides default implementation of the sparse KdTree
 *
 * \see KdTreeDefaultTraits for the default trait interface documentation.
 * \see KdTreeSparseBase for complete API
 */
#ifdef PARSED_WITH_DOXYGEN
/// [KdTreeSparse type definition]
template <typename DataPoint>
struct KdTreeSparse : Ponca::KdTreeSparseBase<KdTreeDefaultTraits<DataPoint>>{};
/// [KdTreeSparse type definition]
#else
template <typename DataPoint>
using KdTreeSparse = KdTreeSparseBase<KdTreeDefaultTraits<DataPoint>>;
#endif

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

    // Queries use a value of -1 for invalid indices
    static_assert(std::is_signed_v<IndexType>, "Index type must be signed");
    static_assert(MAX_DEPTH > 0, "Max depth must be strictly positive");

    // Construction ------------------------------------------------------------
public:
    /// Generate a tree from a custom contained type converted using the specified converter
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam PointConverter Cast/Convert input container type to point container data type
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename PointConverter>
    PONCA_MULTIARCH_HOST inline void build(PointUserContainer&& points, PointConverter c);

    /// Convert a custom point container to the KdTree \ref PointContainer using \ref DataPoint default constructor
    struct DefaultConverter
    {
        template <typename Input>
        PONCA_MULTIARCH_HOST inline void operator()(Input&& i, PointContainer& o)
        {
            using InputContainer = std::remove_reference_t<Input>;
            if constexpr (std::is_same_v<InputContainer, PointContainer> && std::is_copy_assignable_v<typename PointContainer::value_type>)
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
    PONCA_MULTIARCH_HOST inline void build(PointUserContainer&& points)
    {
        build(std::forward<PointUserContainer>(points), DefaultConverter());
    }

    /// Clear tree data
    PONCA_MULTIARCH_HOST inline void clear();

    // Accessors ---------------------------------------------------------------
public:
    PONCA_MULTIARCH [[nodiscard]] inline NodeIndexType nodeCount() const
    {
        return (NodeIndexType)m_nodes_size;
    }

    PONCA_MULTIARCH [[nodiscard]] inline IndexType sampleCount() const
    {
        return (IndexType)m_indices_size;
    }

    PONCA_MULTIARCH [[nodiscard]] inline IndexType pointCount() const
    {
        return (IndexType)m_points_size;
    }

    PONCA_MULTIARCH [[nodiscard]] inline NodeIndexType leafCount() const
    {
        return m_leaf_count;
    }

    PONCA_MULTIARCH [[nodiscard]] inline PointContainer& points()
    {
        return m_points;
    };

    PONCA_MULTIARCH [[nodiscard]] inline const PointContainer& points() const
    {
        return m_points;
    };

    PONCA_MULTIARCH [[nodiscard]] inline const NodeContainer& nodes() const
    {
        return m_nodes;
    }

    PONCA_MULTIARCH [[nodiscard]] inline const IndexContainer& samples() const
    {
        return m_indices;
    }

    // Parameters --------------------------------------------------------------
public:
    /// Read leaf min size
    PONCA_MULTIARCH [[nodiscard]] inline LeafSizeType minCellSize() const
    {
        return m_min_cell_size;
    }

    /// Write leaf min size
    PONCA_MULTIARCH inline void setMinCellSize(LeafSizeType min_cell_size)
    {
        PONCA_DEBUG_ASSERT(min_cell_size > 0);
        m_min_cell_size = min_cell_size;
    }

    // Index mapping -----------------------------------------------------------
public:
    /// Return the point index associated with the specified sample index
    PONCA_MULTIARCH [[nodiscard]] inline IndexType pointFromSample(IndexType sample_index) const
    {
        return m_indices[sample_index];
    }

    /// Return the \ref DataPoint associated with the specified sample index
    /// \note Convenience function, equivalent to
    /// `point_data()[pointFromSample(sample_index)]`
    PONCA_MULTIARCH [[nodiscard]] inline DataPoint& pointDataFromSample(IndexType sample_index)
    {
        return m_points[pointFromSample(sample_index)];
    }
    
    /// Return the \ref DataPoint associated with the specified sample index
    /// \note Convenience function, equivalent to
    /// `point_data()[pointFromSample(sample_index)]`
    PONCA_MULTIARCH [[nodiscard]] inline const DataPoint& pointDataFromSample(IndexType sample_index) const
    {
        return m_points[pointFromSample(sample_index)];
    }

    // Query -------------------------------------------------------------------
public :
    /// \brief Computes a Query object to iterate over the k-nearest neighbors of a point.
    /// The returned object can be reset and reused with the () operator
    /// (using the same argument types as parameters).
    ///
    /// \param point Point from where the query is evaluated
    /// \param k Number of neighbors returned
    /// \return The \ref KdTreeKNearestIndexQuery mutable object to iterate over the search results.
    /// \see KdTreeKNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeKNearestPointQuery<Traits> kNearestNeighbors(const VectorType& point, IndexType k) const
    {
        return KdTreeKNearestPointQuery<Traits>(this, k, point);
    }

    /// \copybrief KdTreeBase::kNearestNeighbors
    /// \param index Index of the point from where the query is evaluated
    /// \param k Number of neighbors returned
    /// \return The \ref KdTreeKNearestIndexQuery mutable object to iterate over the search results.
    /// \see KdTreeKNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeKNearestIndexQuery<Traits> kNearestNeighbors(IndexType index, IndexType k) const
    {
        return KdTreeKNearestIndexQuery<Traits>(this, k, index);
    }

    /// \brief Convenience function that provides an empty k-nearest neighbors Query object.
    ///
    /// The returned object can call for a k-nearest neighbors search using the operator (),
    /// which takes a k and a **position** as parameters.
    ///
    /// Same as `KdTreeBase::kNearestNeighbors (0, VectorType::Zero())`
    /// \return The empty \ref KdTreeKNearestPointQuery mutable object to iterate over the search results.
    /// \see KdTreeKNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeKNearestPointQuery<Traits> kNearestNeighborsQuery() const
    {
        return KdTreeKNearestPointQuery<Traits>(this, 0, VectorType::Zero());
    }

    /// \copybrief KdTreeBase::kNearestNeighborsQuery
    ///
    /// The returned object can call for a k-nearest neighbors search using the operator (),
    /// which takes a k and an **index** as parameters.
    ///
    /// Same as `KdTreeBase::kNearestNeighbors (0, 0)`
    /// \return The empty \ref KdTreeKNearestIndexQuery mutable object to iterate over the search results.
    /// \see KdTreeKNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeKNearestIndexQuery<Traits> kNearestNeighborsIndexQuery() const
    {
        return KdTreeKNearestIndexQuery<Traits>(this, 0, 0);
    }

    /// \brief Computes a Query object that contains the nearest point.
    /// The returned object can be reset and reused with the () operator
    /// (using the same argument types as parameters).
    ///
    /// \param point Point from where the query is evaluated
    /// \return The \ref KdTreeNearestPointQuery mutable object that contains the search result.
    /// \see KdTreeNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeNearestPointQuery<Traits> nearestNeighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<Traits>(this, point);
    }


    /// \copybrief KdTreeBase::nearestNeighbor
    /// \param index Index of the point from where the query is evaluated
    /// \return The \ref KdTreeKNearestIndexQuery mutable object that contains the search result.
    /// \see KdTreeNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeNearestIndexQuery<Traits> nearestNeighbor(IndexType index) const
    {
        return KdTreeNearestIndexQuery<Traits>(this, index);
    }

    /// \brief Convenience function that provides an empty nearest neighbor Query object.
    ///
    /// The returned object can call for a nearest neighbor search using the operator (),
    /// which takes a **position** as parameter.
    ///
    /// Same as `KdTreeBase::nearestNeighbor (VectorType::Zero())`
    ///
    /// \return The empty \ref KdTreeNearestPointQuery mutable object that contains the search result.
    /// \see KdTreeNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeNearestIndexQuery<Traits> nearestNeighborQuery() const
    {
        return KdTreeNearestIndexQuery<Traits>(this, VectorType::Zero());
    }

    /// \copybrief KdTreeBase::nearestNeighborQuery
    ///
    /// The returned object can call for a nearest neighbor search using the operator (),
    /// which takes an **index** as parameter.
    ///
    /// Same as `KdTreeBase::nearestNeighbor (0)`
    ///
    /// \return The \ref KdTreeKNearestIndexQuery mutable object that contains the search result.
    /// \see KdTreeNearestQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeNearestIndexQuery<Traits> nearestNeighborIndexQuery() const
    {
        return KdTreeNearestIndexQuery<Traits>(this, 0);
    }

    /// \brief Computes a Query object to iterate over the neighbors that are inside a given radius.
    /// The returned object can be reset and reused with the () operator
    /// (using the same argument types as parameters).
    ///
    /// \param point Point from where the query is evaluated
    /// \param r Radius around where to search the neighbors
    /// \return The \ref KdTreeRangePointQuery mutable object to iterate over the search results.
    /// \see KdTreeRangeQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeRangePointQuery<Traits> rangeNeighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<Traits>(this, r, point);
    }

    /// \copybrief KdTreeBase::rangeNeighbors
    /// \param index Index of the point from where the query is evaluated
    /// \param r Radius around where to search the neighbors
    /// \return The \ref KdTreeRangeIndexQuery mutable object to iterate over the search results.
    /// \see KdTreeRangeQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeRangeIndexQuery<Traits> rangeNeighbors(IndexType index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<Traits>(this, r, index);
    }

    /// \brief Convenience function that provides an empty range neighbor Query object.
    ///
    /// The returned object can call for a range neighbor search using the operator (),
    /// which takes a **position** as parameter.
    ///
    /// Same as `KdTreeBase::rangeNeighborsQuery (0, VectorType::Zero())`
    ///
    /// \return The empty \ref KdTreeRangePointQuery mutable object to iterate over the search results.
    /// \see KdTreeRangeQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeRangePointQuery<Traits> rangeNeighborsQuery() const
    {
        return KdTreeRangePointQuery<Traits>(this, 0, VectorType::Zero());
    }

    /// \brief KdTreeBase::rangeNeighborsQuery
    ///
    /// The returned object can call for a range neighbor search using the operator (),
    /// which takes an **index** as parameter.
    ///
    /// Same as `KdTreeBase::rangeNeighborsQuery (0, 0)`
    ///
    /// \return The empty \ref KdTreeRangeIndexQuery mutable object to iterate over the search results.
    /// \see KdTreeRangeQueryBase
    PONCA_MULTIARCH [[nodiscard]] KdTreeRangeIndexQuery<Traits> rangeNeighborsIndexQuery() const
    {
        return KdTreeRangeIndexQuery<Traits>(this, 0, 0);
    }

    // Utilities ---------------------------------------------------------------
public:
    PONCA_MULTIARCH_HOST [[nodiscard]] inline bool valid() const;
    PONCA_MULTIARCH_HOST inline void print(std::ostream& os, bool verbose = false) const;

    // Data --------------------------------------------------------------------
protected:
    PointContainer m_points;
    NodeContainer m_nodes;
    IndexContainer m_indices;

    size_t m_points_size{0};
    size_t m_nodes_size{0};
    size_t m_indices_size{0};

    LeafSizeType m_min_cell_size {64}; ///< Minimal number of points per leaf
    NodeIndexType m_leaf_count {0}; ///< Number of leaves in the Kdtree (computed during construction)

    // Internal ----------------------------------------------------------------
public:
    PONCA_MULTIARCH inline KdTreeBase() = default;

    PONCA_MULTIARCH inline KdTreeBase(
        PointContainer points   , NodeContainer nodes    , IndexContainer indices,
        const size_t points_size, const size_t nodes_size, const size_t indices_size
    ) : m_points(points)          , m_nodes(nodes)          , m_indices(indices),
        m_points_size(points_size), m_nodes_size(nodes_size), m_indices_size(indices_size)
    { }
protected:

    /// Generate a tree sampled from a custom contained type converted using a `Converter`
    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam PointConverter Cast/Convert input container type to point container data type
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename PointConverter>
    PONCA_MULTIARCH_HOST inline void buildWithSampling(
        PointUserContainer&& points,
        IndexUserContainer sampling,
        PointConverter c
    );

    /// Generate a tree sampled from a custom contained type converted using a \ref KdTreeBase::DefaultConverter
    /// \tparam PointUserContainer Input points, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Samples used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    PONCA_MULTIARCH_HOST inline void buildWithSampling(PointUserContainer&& points, IndexUserContainer sampling)
    {
        buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling), DefaultConverter());
    }

private:
    PONCA_MULTIARCH_HOST inline void buildRec(std::vector<NodeType>& nodes, NodeIndexType node_id, IndexType start, IndexType end, int level);
    PONCA_MULTIARCH_HOST [[nodiscard]] inline IndexType partition(IndexType start, IndexType end, int dim, Scalar value);
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
    PONCA_MULTIARCH KdTreeDenseBase() = default;

    /// Constructor generating a tree from a custom contained type converted using a \ref Traits::ContainerConverter
    template<typename PointUserContainer>
    PONCA_MULTIARCH_HOST inline explicit KdTreeDenseBase(PointUserContainer&& points)
        : Base()
    {
        Base::build(std::forward<PointUserContainer>(points));
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
    PONCA_MULTIARCH KdTreeSparseBase() = default;

    /// Constructor generating a tree from a custom contained type converted using a \ref Traits::ContainerConverter
    template<typename PointUserContainer>
    PONCA_MULTIARCH_HOST inline explicit KdTreeSparseBase(PointUserContainer&& points)
        : Base()
    {
        this->build(std::forward<PointUserContainer>(points));
    }

    /// Constructor generating a tree sampled from a custom contained type converted using a \ref Traits::ContainerConverter
    /// \tparam PointUserContainer Input points, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Samples used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    PONCA_MULTIARCH_HOST inline KdTreeSparseBase(PointUserContainer&& points, IndexUserContainer sampling)
        : Base()
    {
        this->buildWithSampling(std::forward<PointUserContainer>(points), std::move(sampling));
    }

    using Base::buildWithSampling;
};

#include "./kdTree.hpp"
} // namespace Ponca

template <typename Traits>
PONCA_MULTIARCH_HOST std::ostream& operator<<(std::ostream& os, const Ponca::KdTreeBase<Traits>& kdtree)
{
    kdtree.print(os);
    return os;
}
