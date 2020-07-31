#include <PCA/SpacePartitioning/Grid/Iterator/GridRangeIndexIterator.h>
#include <PCA/SpacePartitioning/Grid/Query/GridRangeIndexQuery.h>

namespace pca {

GridRangeIndexIterator::GridRangeIndexIterator() :
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

GridRangeIndexIterator::GridRangeIndexIterator(GridRangeIndexQuery* query) :
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

GridRangeIndexIterator::GridRangeIndexIterator(GridRangeIndexQuery* query, int index) :
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

bool GridRangeIndexIterator::operator !=(const GridRangeIndexIterator& other) const
{
    return m_index != other.m_index;
}

void GridRangeIndexIterator::operator ++()
{
    m_query->advance(*this);
}

int GridRangeIndexIterator::operator * () const
{
    return m_index;
}

} // namespace pca
