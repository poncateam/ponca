/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <cstddef>

namespace Ponca {

/*!
 *  \brief Input iterator to read the `KdTreeKNearestQueryBase` object.
 *
 *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
 *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator functionalities,
 *  you should copy the index values in an STL-like container.
 *
 *  \note This iterator object can be duplicated with no issues.
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

    PONCA_MULTIARCH inline KdTreeNearestIterator() = default;
    PONCA_MULTIARCH inline KdTreeNearestIterator(Index index) : m_index(index) {}
    PONCA_MULTIARCH virtual inline ~KdTreeNearestIterator() = default;

public:
    /// \brief Inequality operand
    PONCA_MULTIARCH inline bool operator !=(const KdTreeNearestIterator& other) const {
        return m_index != other.m_index;
    }

    /// \brief Equality operand
    PONCA_MULTIARCH inline bool operator ==(const KdTreeNearestIterator& other) const {
        return m_index == other.m_index;
    }

    /// Prefix increment
    PONCA_MULTIARCH inline KdTreeNearestIterator& operator ++() {
        ++m_index;
        return *this;
    }

    /// \brief Postfix increment
    PONCA_MULTIARCH inline KdTreeNearestIterator operator++(int) {
        KdTreeNearestIterator tmp = *this;
        ++m_index;
        return tmp;
    }

    /// \brief Value increment
    PONCA_MULTIARCH inline void operator +=(int i) {m_index += i;}

    /// \brief Dereference operator
    PONCA_MULTIARCH inline reference operator *() const {
        return const_cast<reference>(m_index);
    }

protected:
    Index m_index {-1};
};

} // namespace ponca
