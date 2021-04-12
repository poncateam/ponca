/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template <typename DataPoint>
class KdTreeRangePointIterator
{

public:
    KdTreeRangePointIterator() :
        m_query(nullptr),
        m_index(-1),
        m_start(0),
        m_end(0)
    {
    }

    KdTreeRangePointIterator(KdTreeRangePointQuery<DataPoint>* query) :
        m_query(query),
        m_index(-1),
        m_start(0),
        m_end(0)
    {
    }

    KdTreeRangePointIterator(KdTreeRangePointQuery<DataPoint>* query, int index) :
        m_query(query),
        m_index(index),
        m_start(0),
        m_end(0)
    {
    }

public:
    inline bool operator !=(const KdTreeRangePointIterator<DataPoint>& other) const;
	inline void operator ++(int);
	inline int  operator * () const;

protected:
    KdTreeRangePointQuery<DataPoint>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

#include "./kdTreeRangePointIterator.hpp"
} // namespace ponca
