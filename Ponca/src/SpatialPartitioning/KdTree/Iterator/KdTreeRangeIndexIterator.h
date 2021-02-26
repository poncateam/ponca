#pragma once

#include "../Query/KdTreeRangeIndexQuery.h"

#include "../../query.h"

namespace Ponca {

template<class DataPoint>
class KdTreeRangeIndexIterator
{
protected:
	template<DataPoint> friend class KdTreeRangeIndexQuery;

public:
    KdTreeRangeIndexIterator();
    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery<DataPoint>* query);
    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery<DataPoint>* query, int index);

public:
    bool operator !=(const KdTreeRangeIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    KdTreeRangeIndexQuery<DataPoint>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

} // namespace pdpc
#include "./KdTreeRangeIndexIterator.h"
