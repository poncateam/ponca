/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include <Eigen/Geometry>

#include <cstddef>

#ifdef __has_builtin
#if __has_builtin(__builtin_clz)
#define PONCA_HAS_BUILTIN_CLZ 1
#endif
#endif

#ifndef PONCA_HAS_BUILTIN_CLZ
#define PONCA_HAS_BUILTIN_CLZ 0
#endif

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
    // We're using unsigned indices since we're using bitfields.
    using UIndex = typename std::make_unsigned<NodeIndex>::type;

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
    using Scalar = typename DataPoint::Scalar;

public:
    /*!
     * \brief The type used to store node bounding boxes.
     *
     * Must provide `diagonal()`, and `center()` functions, all returning a
     * `DataPoint::VectorType`.
     */
    using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

    KdTreeDefaultNode() = default;

    /*!*/
    bool is_leaf() const { return m_is_leaf; }

    /*!*/
    void set_is_leaf(bool is_leaf) { m_is_leaf = is_leaf; }

    /*!*/
    void configure_range(Index start, Index size, const AabbType &aabb)
    {
        if (m_is_leaf)
        {
            m_leaf.start = start;
            m_leaf.size = (LeafSize)size;
        }
    }

    /*!*/
    void configure_inner(Scalar split_value, Index first_child_id, Index split_dim)
    {
        if (!m_is_leaf)
        {
            m_inner.split_value = split_value;
            m_inner.first_child_id = first_child_id;
            m_inner.split_dim = split_dim;
        }
    }

    /*!*/
    Index leaf_start() const { return m_leaf.start; }

    /*!*/
    LeafSize leaf_size() const { return m_leaf.size; }

    /*!*/
    Scalar inner_split_value() const { return m_inner.split_value; }
    
    /*!*/
    int inner_split_dim() const { return (int)m_inner.split_dim; }
    
    /*!*/
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
