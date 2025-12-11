/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <iterator>
#include <cstddef>

namespace Ponca {

/*!
 *  \brief Input iterator to read the `KdTreeRangeQuery`.
 *
 *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
 *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator functionalities,
 *  you should copy the values inside a STL-like container.
 *
 *  \warning The increment logic resides in `KdTreeRangeQueryBase::advance(Iterator& it)`
 *  As long as this advance method doesn't update the internal state of the `KdTreeRangeQuery` object,
 *  this iterator can be duplicated without causing issues.
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

    /// \breif Equality operand
    inline bool operator ==(const KdTreeRangeIterator& other) const {
        return m_index == other.m_index;
    }

    /// Prefix increment
    inline KdTreeRangeIterator& operator++() {
        m_query->advance(*this);
        return *this;
    }

    /// \brief Postfix increment
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
