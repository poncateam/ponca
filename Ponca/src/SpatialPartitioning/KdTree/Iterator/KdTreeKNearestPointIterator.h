/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../iterator.h"

namespace Ponca {

template <class DataPoint>
class KdTreeKNearestPointIterator
{
public:
    typedef typename DataPoint::Scalar Scalar; 

    KdTreeKNearestPointIterator() :
        m_iterator()
    {
    }

    KdTreeKNearestPointIterator(typename limited_priority_queue<IndexSquaredDistance<Scalar>>::iterator iterator) :
        m_iterator(iterator)
    {
    }

public:
    bool operator !=(const KdTreeKNearestPointIterator<DataPoint>& other) const;
    void operator ++();
    int  operator * () const;

protected:
    typename limited_priority_queue<IndexSquaredDistance<Scalar>>::iterator m_iterator;
};

#include "./kdTreeKNearestPointIterator.hpp"
} // namespace ponca