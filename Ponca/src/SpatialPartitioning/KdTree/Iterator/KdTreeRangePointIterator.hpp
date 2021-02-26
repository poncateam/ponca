/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "./KdTreeRangePointIterator.h"
#include "../Query/KdTreeRangePointQuery.h"

namespace Ponca {
    /*
template <typename VectorType>
KdTreeRangePointIterator<VectorType>::KdTreeRangePointIterator() :
    m_query(nullptr),
    m_index(-1),
    m_start(0),
    m_end(0)
{
}

template <typename VectorType>
KdTreeRangePointIterator<VectorType>::KdTreeRangePointIterator(KdTreeRangePointQuery<VectorType>* query) :
    m_query(query),
    m_index(-1),
    m_start(0),
    m_end(0)
{
}

template <typename VectorType>
KdTreeRangePointIterator<VectorType>::KdTreeRangePointIterator(KdTreeRangePointQuery<VectorType>* query, int index) :
    m_query(query),
    m_index(index),
    m_start(0),
    m_end(0)
{
}

bool KdTreeRangePointIterator<VectorType>::operator !=(const KdTreeRangePointIterator& other) const
{
    return m_index != other.m_index;
}

void KdTreeRangePointIterator<VectorType>::operator ++()
{
    m_query->advance(*this);
}

int KdTreeRangePointIterator<VectorType>::operator *() const
{
    return m_index;
}*/

} // namespace ponca
