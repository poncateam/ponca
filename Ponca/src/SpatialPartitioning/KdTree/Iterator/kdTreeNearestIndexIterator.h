/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

class KdTreeNearestIndexIterator
{
public:
    KdTreeNearestIndexIterator() :
        m_index(-1)
    {
    }

    KdTreeNearestIndexIterator(int index) :
        m_index(index)
    {
    }

public:
    bool operator !=(const KdTreeNearestIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    int m_index;
};

#include "./kdTreeNearestIndexIterator.hpp"
} // namespace ponca
