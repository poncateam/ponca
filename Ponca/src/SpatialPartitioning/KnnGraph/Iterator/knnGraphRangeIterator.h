/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <cstddef>
#include "../../../Common/Containers/hashset.h"

namespace Ponca
{
    /*
     * - Use `Set = BitSet<MAX_NUMBER_OF_POINTS>` for the fastest search : Provides trivial insertion and search,
     * with O(1) complexity at the expense of memory.
     *
     * - Use `Set = HashSet<Traits::MAX_RANGE_NEIGHBORS_SIZE>` for bigger data set : Best case complexity for insertion
     * and search is O(1) and worst case is O(N) (depends on the given dataset and on the chosen hashing function).
     */
    template <typename Traits, typename Set = HashSet<Traits::MAX_RANGE_NEIGHBORS_SIZE>>
    class KnnGraphRangeQuery;

    /*!
     *  \brief Input iterator to read the `KnnGraphRangeQuery` object.
     *
     *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
     *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator
     * functionalities, you should copy the index values in an STL-like container.
     *
     *  \warning This iterator object should never be duplicated, as it is a proxy that holds a reference to the actual
     * data : The copy of this iterator would still point to the same KnnGraph reference. So, if the increment operator
     * is used on the iterator, the duplicate will also have its state updated. If we then call the increment operator
     * on the duplicate, the result will be an incorrect value.
     *
     *  \see KnnGraphRangeQuery
     */
    template <typename Traits>
    class KnnGraphRangeIterator
    {
    protected:
        friend class KnnGraphRangeQuery<Traits>;
        using Index = typename Traits::IndexType;

    public:
        // Tagged as an input iterator, because the increment logic is shared between iterators of the same queries.
        // Which makes the iterator not valid when duplicated (through postfix increment for example).
        using iterator_category = std::input_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Index;
        using pointer           = Index*;
        using reference         = const Index&;

        PONCA_MULTIARCH inline KnnGraphRangeIterator(KnnGraphRangeQuery<Traits>* query, Index index = Index(-1))
            : m_query(query), m_index(index)
        {
        }

    public:
        /// \brief Inequality operand
        PONCA_MULTIARCH bool operator!=(const KnnGraphRangeIterator& other) const { return m_index != other.m_index; }

        /// \brief Equality operand
        PONCA_MULTIARCH bool operator==(const KnnGraphRangeIterator& other) const { return m_index == other.m_index; }

        /// Prefix increment
        PONCA_MULTIARCH inline KnnGraphRangeIterator& operator++()
        {
            m_query->advance(*this);
            return *this;
        }

        /// \brief Postfix increment
        PONCA_MULTIARCH inline void operator++(value_type) { ++*this; }

        /// \brief Dereference operator
        PONCA_MULTIARCH inline reference operator*() const { return const_cast<reference>(m_index); }

    protected:
        KnnGraphRangeQuery<Traits>* m_query{nullptr};
        value_type m_index{-1};
    };

} // namespace Ponca
