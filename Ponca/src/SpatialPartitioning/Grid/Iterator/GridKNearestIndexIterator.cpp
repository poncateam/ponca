#include <PCA/SpacePartitioning/Grid/Iterator/GridKNearestIndexIterator.h>

namespace pca {

GridKNearestIndexIterator::GridKNearestIndexIterator() :
    m_iterator()
{
}

GridKNearestIndexIterator::GridKNearestIndexIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator) :
    m_iterator(iterator)
{
}

bool GridKNearestIndexIterator::operator !=(const GridKNearestIndexIterator& other) const
{
    return m_iterator != other.m_iterator;
}

void GridKNearestIndexIterator::operator ++()
{
    ++m_iterator;
}

int GridKNearestIndexIterator::operator * () const
{
    return m_iterator->index;
}

} // namespace pca
