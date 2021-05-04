/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

class KdTreeNearestPointIterator
{
public:
    KdTreeNearestPointIterator() :
        m_index(-1)
    {
    }

    KdTreeNearestPointIterator(int index) :
        m_index(index)
    {
    }

public:
    bool operator !=(const KdTreeNearestPointIterator& other) const;
    void operator ++(int);
    inline KdTreeNearestPointIterator& operator++();
    int  operator * () const;

protected:
    int m_index;
};

#include "./kdTreeNearestPointIterator.hpp"
} // namespace ponca
