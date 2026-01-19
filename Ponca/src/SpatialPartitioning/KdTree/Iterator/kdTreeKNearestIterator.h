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
 *  \brief Input iterator to read the `KdTreeKNearestQueryBase` object.
 *
 *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
 *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator functionalities,
 *  you should copy the index values in an STL-like container.
 *
 *  \note This iterator object can be duplicated with no issues.
 *
 *  \see KdTreeKNearestQueryBase
 */
template <typename Index, typename DataPoint>
class KdTreeKNearestIterator
{
public:
    using iterator_category = PONCA_MULTIARCH_CU_STD_NAMESPACE(input_iterator_tag);
    using difference_type   = std::ptrdiff_t;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = const Index&;

    using Scalar   = typename DataPoint::Scalar;
    using Iterator = typename limited_priority_queue<IndexSquaredDistance<Index, Scalar>>::iterator;

    inline KdTreeKNearestIterator() = default;
    inline KdTreeKNearestIterator(const Iterator& iterator) : m_iterator(iterator) {}
    virtual inline ~KdTreeKNearestIterator() = default;

public:
    /// \brief Inequality operand
    inline bool operator !=(const KdTreeKNearestIterator& other) const {
        return m_iterator != other.m_iterator;
    }

    /// \brief Equality operand
    inline bool operator ==(const KdTreeKNearestIterator& other) const {
        return m_iterator == other.m_iterator;
    }

    /// Prefix increment
    inline KdTreeKNearestIterator& operator ++() {++m_iterator; return *this;}

    /// \brief Postfix increment
    inline KdTreeKNearestIterator operator++(int) {
        KdTreeKNearestIterator tmp = *this;
        ++m_iterator;
        return tmp;
    }

    /// \brief Value increment
    inline void operator +=(int i) {m_iterator += i;}

    /// \brief Dereference operator
    inline reference operator *() const {
        return const_cast<reference>(m_iterator->index);
    }

protected:
    Iterator m_iterator;
};
} // namespace ponca
