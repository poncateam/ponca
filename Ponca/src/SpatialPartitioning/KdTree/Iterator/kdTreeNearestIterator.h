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
 *  \brief Input iterator to read the `KdTreeNearestQuery`.
 *
 *  As this is an input iterator, we don't guarantee anything else than reading and incrementing values with it.
 *  If you need to analyse the values with algorithms that relies on forward or more complex iterators,
 *  we suggest copying the values inside a std::vector<Index>.
 *
 *  \note This iterator can be duplicated with no issues.
 *
 *  \see KdTreeNearestQueryBase
 */
template<typename Index>
class KdTreeNearestIterator
{
public:
    using iterator_category = std::input_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = const Index&;

    inline KdTreeNearestIterator() = default;
    inline KdTreeNearestIterator(Index index) : m_index(index) {}
    virtual inline ~KdTreeNearestIterator() = default;

public:
    /// \brief Inequality operand
    inline bool operator !=(const KdTreeNearestIterator& other) const {
        return m_index != other.m_index;
    }

    /// \breif Equality operand
    inline bool operator ==(const KdTreeNearestIterator& other) const {
        return m_index == other.m_index;
    }

    /// Prefix increment
    inline KdTreeNearestIterator& operator ++() {
        ++m_index;
        return *this;
    }

    /// \brief Postfix increment
    inline KdTreeNearestIterator operator++(int) {
        KdTreeNearestIterator tmp = *this;
        ++m_index;
        return tmp;
    }

    /// \brief Value increment
    inline void operator +=(int i) {m_index += i;}

    /// \brief Dereference operator
    inline reference operator *() const {
        return const_cast<reference>(m_index);
    }

protected:
    Index m_index {-1};
};

} // namespace ponca
