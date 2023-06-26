/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include <Eigen/Geometry>

namespace Ponca {
#ifndef PARSED_WITH_DOXYGEN
namespace internal
{
    constexpr int clz(unsigned int value)
    {
#if defined(__has_builtin) && __has_builtin(__builtin_clz)
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

    /*!
     * \brief Calculates bitfield sizes of a default kd-tree inner node at compile-time.
     *
     * \tparam Index The type used to index data points.
     * \tparam DIM The dimension of the data points.
     */
    template <typename Index, int DIM>
    struct KdTreeDefaultInnerNodeBitfieldInfo
    {
        using UIndex = typename std::make_unsigned<Index>::type;

        static_assert(DIM >= 0, "Dim must be positive");

        enum
        {
            /*!
             * \brief The minimum bit width required to store the point dimension.
             *
             * Equal to the index of DIM's most significant bit (starting at 1), e.g.:
             * With DIM = 4,
             *             -------------
             * DIM =    0b | 1 | 0 | 0 |
             *             -------------
             * Bit index =  #3  #2  #1
             *
             * The MSB has an index of 3, so we store the dimension on 3 bits.
             */
            DIM_BITS = sizeof(unsigned int)*8 - internal::clz((unsigned int)DIM),
            UINDEX_BITS = sizeof(UIndex)*8,
        };

        // 1 bit is reserved for the leaf flag
        static_assert(DIM_BITS < UINDEX_BITS - 1,
            "Dim does not fit in the index bitfield of a default inner node");

        enum
        {
            /*!
             * \brief The number of remaining bits that can be used to store indices.
             */
            CHILD_ID_BITS = UINDEX_BITS - (DIM_BITS + 1),
        };
    };
}
#endif

template <typename Index, typename Scalar, int DIM>
struct KdTreeDefaultInnerNode
{
private:
    using BitfieldInfo = internal::KdTreeDefaultInnerNodeBitfieldInfo<Index, DIM>;
    using UIndex       = typename BitfieldInfo::UIndex;

public:
    enum
    {
        /*!
         * \brief The number of bits used to store the point dimension.
         *
         * Points being in higher dimensions will result in less possible points being stored in the kd-tree.
         */
        DIM_BITS = BitfieldInfo::DIM_BITS,

        /*!
         * \brief The number of bits used to store point indices.
         */
        INDEX_BITS = BitfieldInfo::CHILD_ID_BITS,
    };

    Scalar split_value;
    UIndex first_child_id : INDEX_BITS;
    UIndex dim : DIM_BITS;
    UIndex leaf : 1;
};

template <typename Index, typename Size>
struct KdTreeDefaultLeafNode
{
    Index start;
    Size  size;
};

/*!
 * \brief The node type used by default by the kd-tree.
 */
template <typename Index, typename Scalar, int DIM>
struct KdTreeDefaultNode
{
    /*!
     * \brief The type of the `inner` member.
     *
     * Must provide a `split_value` member of type `Scalar`, a `first_child_id` member and a `dim` member.
     *
     * Must also provide `DIM_BITS` and `INDEX_BITS` compile-time constants specifying the number of
     * bits used to store `dim` and `first_child_id` respecitvely. These constants are used to
     * determine the maximum number of points the kd-tree can store.
     */
    using InnerType    = KdTreeDefaultInnerNode<Index, Scalar, DIM>;
    using LeafSizeType = unsigned short;

    /*!
     * \brief The type of the `leaf` member.
     *
     * Must provide a `start` and a `size` member. `size` must be of type `LeafSizeType` and
     * `start` must be large enough to store point indices.
     */
    using LeafType = KdTreeDefaultLeafNode<Index, LeafSizeType>;

    union
    {
        InnerType inner;
        LeafType  leaf;
    };

    bool is_leaf() const { return inner.leaf; }
    void set_is_leaf(bool new_is_leaf) { inner.leaf = new_is_leaf; }
};

/*!
 * \brief The default traits type used by the kd-tree.
 */
template <typename _DataPoint>
struct KdTreeDefaultTraits
{
    /*!
     * \brief The type used to store point data.
     *
     * Must provide `Scalar` and `VectorType` typedefs.
     *
     * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
     * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
     * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
     */
    using DataPoint = _DataPoint;

private:
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;

public:
    enum
    {
        /*!
         * \brief A compile-time constant specifying the maximum depth of the kd-tree.
         */
        MAX_DEPTH = 32,
    };

    /*!
     * \brief The type used to calculate node bounding boxes.
     *
     * Must provide `min()`, `max()`, and `center()` functions, all returning a `VectorType`.
     */
    using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

    // Containers
    using IndexType      = int;
    using PointContainer = std::vector<DataPoint>;
    using IndexContainer = std::vector<IndexType>;

    /*!
     * \brief The container used to store nodes.
     *
     * \note The node type is deduced from the node container type via `NodeContainer::value_type`.
     *
     * \see KdTreeDefaultNode for the node interface documentation.
     */
    using NodeContainer = std::vector<KdTreeDefaultNode<IndexType, Scalar, DataPoint::Dim>>;
};
} // namespace Ponca
