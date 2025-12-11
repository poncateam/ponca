/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template <typename Traits>
class KnnGraphRangeQuery;

/*!
 *  \brief Input iterator to read the `KnnGraphRangeQuery`.
 *
 *  As this is an input iterator, we don't guarantee anything else than reading and incrementing values with it.
 *  If you need to analyse the values with algorithms that relies on forward or more complex iterators,
 *  we suggest copying the values inside a std::vector<Index>.
 *
 *  \warning This iterator should never be duplicated, as it is a proxy that holds a reference to the actual data :
 *  The copy of this iterator would still point to the same KnnGraph reference.
 *  So, if the increment operator is used on the iterator, the duplicate will also have its state updated
 *  in the `KnnGraphRangeQuery`.
 *  If we then call the increment operator on the duplicate, the result will be an incorrect value.
 *
 *  \see KnnGraphRangeQuery::initialize
 */
template <typename Traits>
class KnnGraphRangeIterator
{
protected:
    friend class KnnGraphRangeQuery<Traits>;
    using Index  = typename Traits::IndexType;
public:
    // Tagged as an input iterator, because the increment logic is shared between iterators of the same queries.
    // Which makes the iterator not valid when duplicated (through postfix increment for example).
    using iterator_category = std::input_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = const Index&;

    inline KnnGraphRangeIterator(KnnGraphRangeQuery<Traits>* query, Index index = Index(-1)) : m_query(query), m_index(index) {}

public:
    /// \brief Inequality operand
    bool operator != (const KnnGraphRangeIterator& other) const{
        return m_index != other.m_index;
    }

    /// \brief Equality operand
    bool operator == (const KnnGraphRangeIterator& other) const {
        return m_index == other.m_index;
    }

    /// Prefix increment
    inline KnnGraphRangeIterator& operator ++ (){
        m_query->advance(*this);
        return *this;
    }

    /// \brief Postfix increment
    inline void operator++(value_type) { ++*this; }

    /// \brief Dereference operator
    inline reference operator *() const {
        return const_cast<reference>(m_index);
    }

protected:
    KnnGraphRangeQuery<Traits>* m_query {nullptr};
    value_type m_index {-1};
};

} // namespace Ponca
