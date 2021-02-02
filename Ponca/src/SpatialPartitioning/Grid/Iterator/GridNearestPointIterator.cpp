#include <PCA/SpacePartitioning/Grid/Iterator/GridNearestPointIterator.h>

namespace Ponca {

GridNearestPointIterator::GridNearestPointIterator() :
    m_index(-1)
{
}

GridNearestPointIterator::GridNearestPointIterator(int index) :
    m_index(index)
{
}

bool GridNearestPointIterator::operator !=(const GridNearestPointIterator& other) const
{
    return m_index != other.m_index;
}

void GridNearestPointIterator::operator ++()
{
    ++m_index;
}

int GridNearestPointIterator::operator * () const
{
    return m_index;
}

}   
