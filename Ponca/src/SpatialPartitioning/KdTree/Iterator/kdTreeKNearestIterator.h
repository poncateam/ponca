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
    using Scalar   = typename DataPoint::Scalar;
    using Iterator = typename limited_priority_queue<IndexSquaredDistance<Index, Scalar>>::iterator;

    inline KdTreeKNearestIterator() = default;
    inline KdTreeKNearestIterator(const Iterator& iterator) : m_iterator(iterator) {}
    virtual inline ~KdTreeKNearestIterator() = default;

public:
    inline bool operator !=(const KdTreeKNearestIterator& other) const
    {return m_iterator != other.m_iterator;}
    inline void operator ++() {++m_iterator;}
    inline Index operator * () const {return m_iterator->index;}
    inline void operator +=(int i) {m_iterator += i;}

protected:
    Iterator m_iterator;
};
} // namespace ponca
