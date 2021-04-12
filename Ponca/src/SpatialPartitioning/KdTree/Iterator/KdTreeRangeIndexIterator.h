/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"

namespace Ponca {

template<class DataPoint> class KdTreeRangeIndexQuery;

template<class DataPoint>
class KdTreeRangeIndexIterator
{
protected:
	//friend class KdTreeRangeIndexQuery<DataPoint>;
	friend class KdTreeRangeIndexQuery<DataPoint>;

public:
    KdTreeRangeIndexIterator() :
        m_query(nullptr),
        m_index(-1),
        m_start(0),
        m_end(0)
    {
    }

    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery<DataPoint>* query) :
        m_query(query),
        m_index(-1),
        m_start(0),
        m_end(0)
    {
    }

    KdTreeRangeIndexIterator(KdTreeRangeIndexQuery<DataPoint>* query, int index) :
        m_query(query),
        m_index(index),
        m_start(0),
        m_end(0)
    {
    }

    inline bool operator !=(const KdTreeRangeIndexIterator<DataPoint>& other) const;
    inline void operator ++(int);
    inline int  operator *() const;

protected:
    KdTreeRangeIndexQuery<DataPoint>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

#include "./kdTreeRangeIndexIterator.hpp"
} // namespace ponca
