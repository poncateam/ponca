/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint>
bool KdTreeKNearestPointIterator<DataPoint>::operator !=(const KdTreeKNearestPointIterator<DataPoint>& other) const
{
    return m_iterator != other.m_iterator;
}

template <class DataPoint>
void KdTreeKNearestPointIterator<DataPoint>::operator ++()
{
    ++m_iterator;
}

template <class DataPoint>
int KdTreeKNearestPointIterator<DataPoint>::operator * () const
{
    return m_iterator->index;
}