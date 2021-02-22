#pragma once

#include "../Common/Containers/limitedPriorityQueue.h"
#include "./indexSquaredDistance.h"

namespace Ponca {

struct NearestIterator
{
    inline NearestIterator(int index = -1) : m_index( index ){}

    inline bool operator !=(const NearestIterator& other) const
    { return m_index != other.m_index; }
    inline void operator ++()
    { ++m_index; }
    inline void operator +=(int i)
    { m_index += i; }
    inline int  operator * () const
    { return m_index; }

protected:
    int m_index;
};

//struct KNearestIterator
//{
//    inline KNearestIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator
//                            = limited_priority_queue<IndexSquaredDistance>::iterator() )
//        : m_iterator(iterator) {}
//
//    inline bool operator !=(const KNearestIterator& other) const
//    { return m_iterator != other.m_iterator; }
//    inline void operator ++()
//    { ++m_iterator; }
//    inline void operator +=(int i)
//    { m_iterator += i; }
//    inline int  operator * () const
//    { return m_iterator->index; }
//
//
//protected:
//    limited_priority_queue<IndexSquaredDistance>::iterator m_iterator;
//};

template<typename _RangeQueryType>
struct RangeIterator
{
    using RangeQueryType = _RangeQueryType;

    inline RangeIterator(RangeQueryType* query = nullptr)
        : m_query( query ) {}
    inline RangeIterator(RangeQueryType* query, int index)
        : m_index( index ), m_query( query ) {}

    inline bool operator !=(const RangeIterator& other) const
    { return m_index != other.m_index; }
    inline int  operator * () const
    { return m_index; }
    inline void operator++ () const { m_query->advance(*this); }

protected:
    int m_index {-1};
    int m_start { 0};
    int m_end   { 0};
    RangeQueryType * m_query;
};

} // namespace Ponca