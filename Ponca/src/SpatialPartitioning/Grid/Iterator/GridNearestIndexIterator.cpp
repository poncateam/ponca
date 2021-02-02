#include <PCA/SpacePartitioning/Grid/Iterator/GridNearestIndexIterator.h>

namespace Ponca {

GridNearestIndexIterator::GridNearestIndexIterator() :
    m_index(-1)
{
}

GridNearestIndexIterator::GridNearestIndexIterator(int index) :
    m_index(index)
{
}

bool GridNearestIndexIterator::operator !=(const GridNearestIndexIterator& other) const
{
    return m_index != other.m_index;
}

void GridNearestIndexIterator::operator ++()
{
    ++m_index;
}

int GridNearestIndexIterator::operator * () const
{
    return m_index;
}

}   
