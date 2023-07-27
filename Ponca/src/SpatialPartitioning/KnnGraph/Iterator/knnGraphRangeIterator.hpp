/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "KnnGraphRangeIterator.h"
#include "../Query/knnGraphRangeQuery.h"

namespace Ponca {

KnnGraphRangeIterator::KnnGraphRangeIterator() :
    m_query(nullptr),
    m_index(-1)
{
}

KnnGraphRangeIterator::KnnGraphRangeIterator(KnnGraphRangeQuery* query) :
    m_query(query),
    m_index(-1)
{
}

KnnGraphRangeIterator::KnnGraphRangeIterator(KnnGraphRangeQuery* query, int index) :
    m_query(query),
    m_index(index)
{
}

bool KnnGraphRangeIterator::operator != (const KnnGraphRangeIterator& other) const
{
    return m_index != other.m_index;
}

void KnnGraphRangeIterator::operator ++ ()
{
    m_query->advance(*this);
}

int KnnGraphRangeIterator::operator *  () const
{
    return m_index;
}


} // namespace Ponca
