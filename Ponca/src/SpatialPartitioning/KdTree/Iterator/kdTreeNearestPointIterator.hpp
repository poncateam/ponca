/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

bool KdTreeNearestPointIterator::operator !=(const KdTreeNearestPointIterator& other) const
{
    return m_index != other.m_index;
}

void KdTreeNearestPointIterator::operator ++()
{
    ++m_index;
}

int KdTreeNearestPointIterator::operator * () const
{
    return m_index;
}

