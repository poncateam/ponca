/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

bool KdTreeNearestIndexIterator::operator !=(const KdTreeNearestIndexIterator& other) const
{
    return m_index != other.m_index;
}

void KdTreeNearestIndexIterator::operator ++(int)
{
    ++m_index;
}

KdTreeNearestIndexIterator& KdTreeNearestIndexIterator::operator ++()
{
    ++m_index;
    return *this;
}

int KdTreeNearestIndexIterator::operator * () const
{
    return m_index;
}
