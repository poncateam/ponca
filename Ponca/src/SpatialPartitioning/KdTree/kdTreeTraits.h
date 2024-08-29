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

    Scalar split_value{0};
    UIndex first_child_id : INDEX_BITS;
    UIndex split_dim : DIM_BITS;
};

template <typename Index, typename Size>
struct KdTreeDefaultLeafNode
{
    Index start{0};
    Size size{0};
};

/*!
 * \brief The node type used by default by the kd-tree.
 *
 * It is possible to modify the Inner and Leaf node types by inheritance. For instance, to add a Bounding box to inner
 * nodes, define a custom inner node type:
 *
 * \snippet ponca_customize_kdtree.cpp CustomInnerNodeDefinition
 *
 * Define a custom node type to use it, and expose custom data (inner/leaf node are not exposed directly):
 *
 * \snippet ponca_customize_kdtree.cpp CustomNodeDefinition
 *
 * To use in the KdTree, define a type using the custom node:
 *
 * \snippet ponca_customize_kdtree.cpp KdTreeTypeWithCustomNode
 *
 * The added attribute can be accessed
 *
 * \snippet ponca_customize_kdtree.cpp ReadCustomProperties
 *
 */
template <typename Index, typename NodeIndex, typename DataPoint,
          typename LeafSize = Index,
          typename _InnerNodeType = KdTreeDefaultInnerNode<NodeIndex, typename DataPoint::Scalar, DataPoint::Dim>,
          typename _LeafNodeType = KdTreeDefaultLeafNode<Index, LeafSize>>
class KdTreeCustomizableNode
{
private:
    using Scalar    = typename DataPoint::Scalar;
    using InnerType = _InnerNodeType;
    using LeafType  = _LeafNodeType;

public:
    enum : std::size_t
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
    
    [[nodiscard]] bool is_leaf() const { return m_is_leaf; }
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
            data.m_leaf.start = start;
            data.m_leaf.size = (LeafSize)size;
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
            data.m_inner.split_value = split_value;
            data.m_inner.first_child_id = first_child_id;
            data.m_inner.split_dim = split_dim;
        }
    }

    /*!
     * \brief The start index of the range of the leaf node in the sample
     * index array.
     */
    [[nodiscard]] Index leaf_start() const { return data.m_leaf.start; }

    /*!
     * \brief The size of the range of the leaf node in the sample index array.
     */
    [[nodiscard]] LeafSize leaf_size() const { return data.m_leaf.size; }

    /*!
     * \brief The position of the AABB split of the inner node.
     */
    [[nodiscard]] Scalar inner_split_value() const { return data.m_inner.split_value; }
    
    /*!
     * \brief Which axis the split of the AABB of the inner node was done on.
     */
    [[nodiscard]] int inner_split_dim() const { return (int)data.m_inner.split_dim; }
    
    /*!
     * \brief The index of the first child of the node in the node array of the
     * kd-tree.
     *
     * \note The second child is stored directly after the first in the array
     * (i.e. `first_child_id + 1`).
     */
    [[nodiscard]] Index inner_first_child_id() const { return (Index)data.m_inner.first_child_id; }

protected:
    [[nodiscard]] inline LeafType& getAsLeaf() { return data.m_leaf; }
    [[nodiscard]] inline InnerType& getAsInner() { return data.m_inner; }
    [[nodiscard]] inline const LeafType& getAsLeaf() const { return data.m_leaf; }
    [[nodiscard]] inline const InnerType& getAsInner() const { return data.m_inner; }

private:
    bool m_is_leaf{true};
    union Data
    {
        // We need an explicit constructor here, see https://stackoverflow.com/a/70428826
        constexpr Data() : m_leaf() {}
        // Needed to satisfy MoveInsertable requirement https://en.cppreference.com/w/cpp/named_req/MoveInsertable
        constexpr Data(const Data&d) : m_leaf(d.m_leaf) {}

        ~Data() {}
        LeafType m_leaf;
        InnerType m_inner;
    };
    Data data;
};

template <typename Index, typename NodeIndex, typename DataPoint,
        typename LeafSize = Index>
struct KdTreeDefaultNode : public KdTreeCustomizableNode<Index, NodeIndex, DataPoint, LeafSize,
        KdTreeDefaultInnerNode<NodeIndex, typename DataPoint::Scalar, DataPoint::Dim>,
        KdTreeDefaultLeafNode<Index, LeafSize>> {
    using Base =  KdTreeCustomizableNode<Index, NodeIndex, DataPoint, LeafSize,
            KdTreeDefaultInnerNode<NodeIndex, typename DataPoint::Scalar, DataPoint::Dim>,
            KdTreeDefaultLeafNode<Index, LeafSize>>;
};

/*!
 * \brief The default traits type used by the kd-tree.
 *
 * \see KdTreeCustomizableNode Helper class to modify Inner/Leaf nodes without redefining a Trait class
 *
 * \tparam _NodeType Type used to store nodes, set by default to #KdTreeDefaultNode
 */
template <typename _DataPoint,
        template <typename /*Index*/,
                  typename /*NodeIndex*/,
                  typename /*DataPoint*/,
                  typename /*LeafSize*/> typename _NodeType = KdTreeDefaultNode>
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
    using NodeType      = _NodeType<IndexType, NodeIndexType, DataPoint, LeafSizeType>;
    using NodeContainer = std::vector<NodeType>;
};
} // namespace Ponca
