/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

//#include "../Query/KdTreeRangePointQuery.h"

//#include "../../query.h"


template<class DataPoint> class KdTreeRangePointQuery;

namespace Ponca {

template <typename DataPoint>
class KdTreeRangePointIterator
{
public:
    /*! \brief Scalar type inherited from DataPoint */
    typedef typename DataPoint::Scalar     Scalar;
    /*! \brief Vector type inherited from DataPoint */
    typedef typename DataPoint::VectorType VectorType;
protected:
    friend class KdTreeRangePointQuery<DataPoint>;

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
    bool operator !=(const KdTreeRangePointIterator<DataPoint>& other) const;
    void operator ++();
    int  operator * () const;

protected:
    KdTreeRangePointQuery<DataPoint>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

#include "./KdTreeRangePointIterator.hpp"
} // namespace ponca
