/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../iterator.h"

namespace Ponca {

class KdTreeKNearestPointIterator
{
public:
    KdTreeKNearestPointIterator();
    KdTreeKNearestPointIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator);

public:
    bool operator !=(const KdTreeKNearestPointIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    limited_priority_queue<IndexSquaredDistance>::iterator m_iterator;
};

} // namespace ponca
