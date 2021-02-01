#include <PCA/SpacePartitioning/KNNGraph/Iterator/KNNGraphRangeIterator.h>
#include <PCA/SpacePartitioning/KNNGraph/Query/KNNGraphRangeQuery.h>

namespace Ponca {

KNNGraphRangeIterator::KNNGraphRangeIterator() :
    m_query(nullptr),
    m_index(-1)
{
}

KNNGraphRangeIterator::KNNGraphRangeIterator(KNNGraphRangeQuery* query) :
    m_query(query),
    m_index(-1)
{
}

KNNGraphRangeIterator::KNNGraphRangeIterator(KNNGraphRangeQuery* query, int index) :
    m_query(query),
    m_index(index)
{
}

bool KNNGraphRangeIterator::operator != (const KNNGraphRangeIterator& other) const
{
    return m_index != other.m_index;
}

void KNNGraphRangeIterator::operator ++ ()
{
    m_query->advance(*this);
}

int KNNGraphRangeIterator::operator *  () const
{
    return m_index;
}


} // namespace pca
