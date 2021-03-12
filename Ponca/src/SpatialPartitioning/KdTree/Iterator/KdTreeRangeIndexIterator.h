/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

//#include "../Query/KdTreeRangeIndexQuery.h"

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
    inline void operator ++();
    inline int  operator * () const;

	/*
	inline bool KdTreeRangeIndexIterator<DataPoint>::operator !=(const KdTreeRangeIndexIterator<DataPoint>& other) const
	{
		return m_index != other.m_index;
	}

	inline void KdTreeRangeIndexIterator<DataPoint>::operator ++()
	{
		m_query->advance(*this);
	}

	inline int KdTreeRangeIndexIterator<DataPoint>::operator *() const
	{
		return m_index;
	}*/

protected:
    KdTreeRangeIndexQuery<DataPoint>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

#include "./KdTreeRangeIndexIterator.hpp"
} // namespace ponca
