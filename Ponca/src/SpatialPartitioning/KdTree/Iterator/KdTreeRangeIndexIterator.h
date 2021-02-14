#pragma once

#include "../Query/KdTreeRangeIndexQuery.h"

#include "../../query.h"

namespace Ponca {

class KdTreeRangeIndexQuery;

class KdTreeRangeIndexIterator
{
protected:
    friend class KdTreeRangeIndexQuery;

public:
    KdTreeRangeIndexIterator();
    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery* query);
    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery* query, int index);

public:
    bool operator !=(const KdTreeRangeIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    KdTreeRangeIndexQuery* m_query;
    int m_index;
    int m_start;
    int m_end;
};

} // namespace pdpc
#include "./KdTreeRangeIndexIterator.h"
