/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template <typename Index, typename DataPoint>
class KdTreeKNearestIterator
{
public:
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = Index&;

    using Scalar   = typename DataPoint::Scalar;
    using Iterator = typename limited_priority_queue<IndexSquaredDistance<Index, Scalar>>::iterator;

    inline KdTreeKNearestIterator() = default;
    inline KdTreeKNearestIterator(const Iterator& iterator) : m_iterator(iterator) {}
    virtual inline ~KdTreeKNearestIterator() = default;

public:
    inline bool operator !=(const KdTreeKNearestIterator& other) const
    {return m_iterator != other.m_iterator;}
    inline void operator ++() {++m_iterator;}
    inline reference operator *() const {
        return const_cast<reference>(m_iterator->index);
    }
    inline void operator +=(int i) {m_iterator += i;}

protected:
    Iterator m_iterator;
};
} // namespace ponca
