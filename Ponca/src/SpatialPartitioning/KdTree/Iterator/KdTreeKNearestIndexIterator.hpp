/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint>
bool KdTreeKNearestIndexIterator<DataPoint>::operator !=(const KdTreeKNearestIndexIterator<DataPoint>& other) const
{
    return m_iterator != other.m_iterator;
}

template <class DataPoint>
void KdTreeKNearestIndexIterator<DataPoint>::operator ++()
{
    ++m_iterator;
}

template <class DataPoint>
int KdTreeKNearestIndexIterator<DataPoint>::operator * () const
{
    return m_iterator->index;
}

template <class DataPoint>
void KdTreeKNearestIndexIterator<DataPoint>::operator +=(int i)
{
    m_iterator += i;
}
