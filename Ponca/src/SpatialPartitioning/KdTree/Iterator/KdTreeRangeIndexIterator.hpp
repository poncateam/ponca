
#include "./KdTreeRangeIndexIterator.h"
#include "../Query/KdTreeRangePointQuery.h"

namespace Ponca {
/*
KdTreeRangeIndexIterator::KdTreeRangeIndexIterator() :
    m_query(nullptr),
    m_index(-1),
    m_start(0),
    m_end(0)
{
}

KdTreeRangeIndexIterator::KdTreeRangeIndexIterator(KdTreeRangeIndexQuery* query) :
    m_query(query),
    m_index(-1),
    m_start(0),
    m_end(0)
{
}

KdTreeRangeIndexIterator::KdTreeRangeIndexIterator(KdTreeRangeIndexQuery* query, int index) :
    m_query(query),
    m_index(index),
    m_start(0),
    m_end(0)
{
}

bool KdTreeRangeIndexIterator::operator !=(const KdTreeRangeIndexIterator& other) const
{
    return m_index != other.m_index;
}

void KdTreeRangeIndexIterator::operator ++()
{
    m_query->advance(*this);
}

int KdTreeRangeIndexIterator::operator * () const
{
    return m_index;
}*/

} // namespace pdpc
