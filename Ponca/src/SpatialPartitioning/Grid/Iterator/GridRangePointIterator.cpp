#include <PCA/SpacePartitioning/Grid/Iterator/GridRangePointIterator.h>
#include <PCA/SpacePartitioning/Grid/Query/GridRangePointQuery.h>

namespace pca {

GridRangePointIterator::GridRangePointIterator() :
    m_query(nullptr),
    m_index(-1),
    m_i_start(-1),
    m_j_start(-1),
    m_k_start(-1),
    m_i_end(-1),
    m_j_end(-1),
    m_k_end(-1),
    m_i(-1),
    m_j(-1),
    m_k(-1),
    m_it(-1),
    m_idx_cell(-1)
{
}

GridRangePointIterator::GridRangePointIterator(GridRangePointQuery* query) :
    m_query(query),
    m_index(-1),
    m_i_start(-1),
    m_j_start(-1),
    m_k_start(-1),
    m_i_end(-1),
    m_j_end(-1),
    m_k_end(-1),
    m_i(-1),
    m_j(-1),
    m_k(-1),
    m_it(-1),
    m_idx_cell(-1)
{
}

GridRangePointIterator::GridRangePointIterator(GridRangePointQuery* query, int index) :
    m_query(query),
    m_index(index),
    m_i_start(-1),
    m_j_start(-1),
    m_k_start(-1),
    m_i_end(-1),
    m_j_end(-1),
    m_k_end(-1),
    m_i(-1),
    m_j(-1),
    m_k(-1),
    m_it(-1),
    m_idx_cell(-1)
{
}

bool GridRangePointIterator::operator !=(const GridRangePointIterator& other) const
{
    return m_index != other.m_index;
}

void GridRangePointIterator::operator ++()
{
    m_query->advance(*this);
}

int GridRangePointIterator::operator * () const
{
    return m_index;
}

} // namespace pca
