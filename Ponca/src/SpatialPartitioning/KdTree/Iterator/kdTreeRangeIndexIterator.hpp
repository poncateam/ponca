/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint>
bool KdTreeRangeIndexIterator<DataPoint>::operator !=(const KdTreeRangeIndexIterator<DataPoint>& other) const
{
    return m_index != other.m_index;
}

template<class DataPoint>
void KdTreeRangeIndexIterator<DataPoint>::operator ++(int)
{
    m_query->advance(*this);
}

template<class DataPoint>
int KdTreeRangeIndexIterator<DataPoint>::operator *() const
{
    return m_index;
}



