/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include "../../Common/Macro.h"

#include <cstddef>

#include <Eigen/Geometry>

namespace Ponca {
#ifndef PARSED_WITH_DOXYGEN
namespace internal
{
    constexpr int clz(unsigned int value)
    {
#if PONCA_HAS_BUILTIN_CLZ
        return __builtin_clz(value);
#else
        if (value == 0)
        {
            return 0;
        }

        unsigned int msb_mask = 1 << (sizeof(unsigned int)*8 - 1);
        int count = 0;
        for (; (value & msb_mask) == 0; value <<= 1, ++count)
        {
        }
        return count;
#endif
    }
}
#endif

template <typename NodeIndex, typename Scalar, int DIM>
struct KdTreeDefaultInnerNode
{
private:
    enum
    {
        // The minimum bit width required to store the split dimension.
        // 
        // Equal to the index of DIM's most significant bit (starting at 1), e.g.:
        // With DIM = 4,
        //             -------------
        // DIM =    0b | 1 | 0 | 0 |
        //             -------------
        // Bit index =  #3  #2  #1
        // 
        // The MSB has an index of 3, so we store the dimension on 3 bits.
        DIM_BITS = sizeof(unsigned int)*8 - internal::clz((unsigned int)DIM),
    };

    // The node stores bitfields as unsigned indices.
    using UIndex = typename std::make_unsigned<NodeIndex>::type;

public:
    enum
    {
        /*!
         * \brief The bit width used to store the first child index.
         */
        INDEX_BITS = sizeof(UIndex)*8 - DIM_BITS,
    };

    Scalar split_value;
    UIndex first_child_id : INDEX_BITS;
    UIndex split_dim : DIM_BITS;
};

template <typename Index, typename Size>
struct KdTreeDefaultLeafNode
{
    Index start;
    Size size;
};

/*!
 * \brief The node type used by default by the kd-tree.
 */
template <typename Index, typename NodeIndex, typename DataPoint,
          typename LeafSize = Index>
class KdTreeDefaultNode
{
private:
    using Scalar    = typename DataPoint::Scalar;
    using LeafType  = KdTreeDefaultLeafNode<Index, LeafSize>;
    using InnerType = KdTreeDefaultInnerNode<NodeIndex, Scalar, DataPoint::Dim>;

public:
    enum
    {
        /*!
         * \brief The maximum number of nodes that a kd-tree can have when using
         * this node type.
         */
        MAX_COUNT = std::size_t(2) << InnerType::INDEX_BITS,
    };

    /*!
     * \brief The type used to store node bounding boxes.
     *
     * Must provide `diagonal()`, and `center()` functions, all returning a
     * `DataPoint::VectorType`.
     */
    using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

    KdTreeDefaultNode() = default;
    
    bool is_leaf() const { return m_is_leaf; }
    void set_is_leaf(bool is_leaf) { m_is_leaf = is_leaf; }

    /*!
     * \brief Configures the range of the node in the sample index array of the
     * kd-tree.
     *
     * \see the leaf node accessors for a more detailed explanation of each
     * argument.
     *
     * \note The AABB is not required by the implementation, so nodes don't
     * have to make it available.
     *
     * Called after \ref set_is_leaf during kd-tree construction.
     */
    void configure_range(Index start, Index size, const AabbType &aabb)
    {
        if (m_is_leaf)
        {
            m_leaf.start = start;
            m_leaf.size = (LeafSize)size;
        }
    }

    /*!
     * \brief Configures the inner node information.
     *
     * \see the inner node accessors for a more detailed explanation of each
     * argument.
     *
     * Called after \ref set_is_leaf and \ref configure_range during kd-tree
     * construction.
     */
    void configure_inner(Scalar split_value, Index first_child_id, Index split_dim)
    {
        if (!m_is_leaf)
        {
            m_inner.split_value = split_value;
            m_inner.first_child_id = first_child_id;
            m_inner.split_dim = split_dim;
        }
    }

    /*!
     * \brief The start index of the range of the leaf node in the sample
     * index array.
     */
    Index leaf_start() const { return m_leaf.start; }

    /*!
     * \brief The size of the range of the leaf node in the sample index array.
     */
    LeafSize leaf_size() const { return m_leaf.size; }

    /*!
     * \brief The position of the AABB split of the inner node.
     */
    Scalar inner_split_value() const { return m_inner.split_value; }
    
    /*!
     * \brief Which axis the split of the AABB of the inner node was done on.
     */
    int inner_split_dim() const { return (int)m_inner.split_dim; }
    
    /*!
     * \brief The index of the first child of the node in the node array of the
     * kd-tree.
     *
     * \note The second child is stored directly after the first in the array
     * (i.e. `first_child_id + 1`).
     */
    Index inner_first_child_id() const { return (Index)m_inner.first_child_id; }

private:
    bool m_is_leaf;
    union
    {
        KdTreeDefaultLeafNode<Index, LeafSize> m_leaf;
        KdTreeDefaultInnerNode<NodeIndex, Scalar, DataPoint::Dim> m_inner;
    };
};

/*!
 * \brief The default traits type used by the kd-tree.
 */
template <typename _DataPoint>
struct KdTreeDefaultTraits
{
    enum
    {
        /*!
         * \brief A compile-time constant specifying the maximum depth of the kd-tree.
         */
        MAX_DEPTH = 32,
    };

    /*!
     * \brief The type used to store point data.
     *
     * Must provide `Scalar` and `VectorType` typedefs.
     *
     * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
     * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
     * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
     */
    using DataPoint    = _DataPoint;
    using IndexType    = int;
    using LeafSizeType = unsigned short;

    // Containers
    using PointContainer = std::vector<DataPoint>;
    using IndexContainer = std::vector<IndexType>;

    // Nodes
    using NodeIndexType = std::size_t;
    using NodeType      = KdTreeDefaultNode<IndexType, NodeIndexType, DataPoint, LeafSizeType>;
    using NodeContainer = std::vector<NodeType>;
};
} // namespace Ponca
