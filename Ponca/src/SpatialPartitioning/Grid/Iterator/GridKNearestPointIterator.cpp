#include <PCA/SpacePartitioning/Grid/Iterator/GridKNearestPointIterator.h>

namespace Ponca {

GridKNearestPointIterator::GridKNearestPointIterator() :
    m_iterator()
{
}

GridKNearestPointIterator::GridKNearestPointIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator) :
    m_iterator(iterator)
{
}

bool GridKNearestPointIterator::operator !=(const GridKNearestPointIterator& other) const
{
    return m_iterator != other.m_iterator;
}

void GridKNearestPointIterator::operator ++()
{
    ++m_iterator;
}

int GridKNearestPointIterator::operator * () const
{
    return m_iterator->index;
}

}   
