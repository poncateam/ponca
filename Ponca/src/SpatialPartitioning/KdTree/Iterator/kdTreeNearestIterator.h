/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

class KdTreeNearestIterator
{
public:
    inline KdTreeNearestIterator() = default;
    inline KdTreeNearestIterator(int index) : m_index(index) {}
    virtual inline ~KdTreeNearestIterator() = default;

public:
    inline bool operator !=(const KdTreeNearestIterator& other) const
    {return m_index != other.m_index;}
    inline void operator ++(int) {++m_index;}
    inline KdTreeNearestIterator& operator ++() {++m_index; return *this;}
    inline int  operator * () const {return m_index;}

protected:
    int m_index {-1};
};

} // namespace ponca
