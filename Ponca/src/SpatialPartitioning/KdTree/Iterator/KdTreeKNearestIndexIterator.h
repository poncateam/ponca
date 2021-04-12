/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template <class DataPoint>
class KdTreeKNearestIndexIterator
{
public:
    typedef typename DataPoint::Scalar Scalar;
    typedef typename limited_priority_queue<IndexSquaredDistance<Scalar>>::iterator Iterator;

    KdTreeKNearestIndexIterator() :
        m_iterator()
    {
    }

    KdTreeKNearestIndexIterator(Iterator iterator) :
        m_iterator(iterator)
    {
    }

public:
    bool operator !=(const KdTreeKNearestIndexIterator<DataPoint>& other) const;
    void operator ++();
    int  operator * () const;
    void operator +=(int i);

protected:
    Iterator m_iterator;
};

#include "./kdTreeKNearestIndexIterator.hpp"
} // namespace ponca
