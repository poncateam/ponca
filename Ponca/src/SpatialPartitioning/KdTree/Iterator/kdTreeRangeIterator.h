/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <cstddef>

namespace Ponca {

/*!
 *  \brief Input iterator to read the `KdTreeRangeQueryBase` object.
 *
 *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
 *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator functionalities,
 *  you should copy the index values in an STL-like container.
 *
 *  \note The increment logic resides in `KdTreeRangeQueryBase::advance(Iterator& it)`
 *  Since the `KdTreeRangeQueryBase::advance` method doesn't update the internal state of the
 *  `KdTreeRangeQuery` object, this iterator object can be duplicated without causing issues.
 *  If a copy of this iterator is made  (e.g., passed by copy to a function),
 *  incrementing one iterator won't update the state of the other.
 *
 *  \see KdTreeRangeQueryBase
 */
template<typename Index, typename DataPoint, typename QueryT_>
class KdTreeRangeIterator
{
protected:
    friend QueryT_;

public:
    using iterator_category = std::input_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = const Index&;

    using Scalar     = typename DataPoint::Scalar;
    using QueryType  = QueryT_;

    inline KdTreeRangeIterator() = default;
    inline KdTreeRangeIterator(QueryType* query, Index index = -1) :
        m_query(query), m_index(index), m_start(0), m_end(0) {}

    /// \brief Inequality operand
    inline bool operator !=(const KdTreeRangeIterator& other) const {
        return m_index != other.m_index;
    }

    /// \brief Equality operand
    inline bool operator ==(const KdTreeRangeIterator& other) const {
        return m_index == other.m_index;
    }

    /// \brief Prefix increment
    /// \see KdTreeRangeQueryBase::advance(Iterator& it) for the iteration logic
    inline KdTreeRangeIterator& operator++() {
        m_query->advance(*this);
        return *this;
    }

    /// \brief Postfix increment
    /// \see KdTreeRangeQueryBase::advance(Iterator& it) for the iteration logic
    inline KdTreeRangeIterator operator++(int) {
        KdTreeRangeIterator tmp = *this;
        m_query->advance(*this);
        return tmp;
    }

    /// \brief Dereference operator
    inline reference operator *() const {
       return const_cast<reference>(m_index);
    }

protected:
    QueryType* m_query {nullptr};
    Index m_index {-1};
    Index m_start {0};
    Index m_end {0};
};
} // namespace ponca
