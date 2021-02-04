/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Query/KdTreeRangePointQuery.h"

#include "../../query.h"

namespace Ponca {

template <typename VectorType>
class KdTreeRangePointIterator
{
protected:
    template <typename _VectorType>
    friend class KdTreeRangePointQuery;

public:
    KdTreeRangePointIterator();
    KdTreeRangePointIterator(KdTreeRangePointQuery<VectorType>* query);
    KdTreeRangePointIterator(KdTreeRangePointQuery<VectorType>* query, int index);

public:
    bool operator !=(const KdTreeRangePointIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    KdTreeRangePointQuery<VectorType>* m_query;
    int m_index;
    int m_start;
    int m_end;
};

} // namespace ponca
#include "./KdTreeRangePointIterator.hpp"
