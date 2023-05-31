/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template<typename Index, typename DataPoint, typename QueryT_>
class KdTreeRangeIterator
{
protected:
    friend QueryT_;

public:
    using Scalar    = typename DataPoint::Scalar;
    using QueryType = QueryT_;

    inline KdTreeRangeIterator() = default;
    inline KdTreeRangeIterator(QueryType* query, Index index = -1) :
        m_query(query), m_index(index), m_start(0), m_end(0) {}

    inline bool operator !=(const KdTreeRangeIterator& other) const
    {return m_index != other.m_index;}
    inline void operator ++(int) {m_query->advance(*this);}
    inline KdTreeRangeIterator& operator++() {m_query->advance(*this); return *this;}
    inline Index operator *() const {return m_index;}

protected:
    QueryType* m_query {nullptr};
    Index m_index {-1};
    Index m_start {0};
    Index m_end {0};
};
} // namespace ponca
