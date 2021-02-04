/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "./KdTreeKNearestPointIterator.h"

namespace Ponca {

KdTreeKNearestPointIterator::KdTreeKNearestPointIterator() :
    m_iterator()
{
}

KdTreeKNearestPointIterator::KdTreeKNearestPointIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator) :
    m_iterator(iterator)
{
}

bool KdTreeKNearestPointIterator::operator !=(const KdTreeKNearestPointIterator& other) const
{
    return m_iterator != other.m_iterator;
}

void KdTreeKNearestPointIterator::operator ++()
{
    ++m_iterator;
}

int KdTreeKNearestPointIterator::operator * () const
{
    return m_iterator->index;
}

} // namespace ponca
