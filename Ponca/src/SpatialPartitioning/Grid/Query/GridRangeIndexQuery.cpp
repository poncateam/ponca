#include <PCA/SpacePartitioning/Grid/Query/GridRangeIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Grid.h>

namespace pca {

GridRangeIndexQuery::GridRangeIndexQuery() :
    GridQuery(),
    RangeIndexQuery()
{
}

GridRangeIndexQuery::GridRangeIndexQuery(const Grid* grid) :
    GridQuery(grid),
    RangeIndexQuery()
{
}

GridRangeIndexQuery::GridRangeIndexQuery(const Grid* grid, Scalar radius) :
    GridQuery(grid),
    RangeIndexQuery(radius)
{
}

GridRangeIndexQuery::GridRangeIndexQuery(const Grid* grid, Scalar radius, int index) :
    GridQuery(grid),
    RangeIndexQuery(radius, index)
{
}

GridRangeIndexIterator GridRangeIndexQuery::begin()
{
    GridRangeIndexIterator it(this);
    this->initialize(it);
    this->advance(it);
    return it;
}

GridRangeIndexIterator GridRangeIndexQuery::end()
{
    return GridRangeIndexIterator(this, m_grid->size());
}

void GridRangeIndexQuery::initialize(GridRangeIndexIterator& it)
{
    const auto& point = m_grid->point_data()[m_index];
    const auto radius = this->radius();
    const auto idx3D = m_grid->cell_index_3D(point);

    int cell_count = std::ceil(radius / m_grid->cell_size());

    int imin = std::max(0, idx3D.i - cell_count);
    int jmin = std::max(0, idx3D.j - cell_count);
    int kmin = std::max(0, idx3D.k - cell_count);
    int imax = std::min(m_grid->cell_count_x()-1, idx3D.i + cell_count);
    int jmax = std::min(m_grid->cell_count_y()-1, idx3D.j + cell_count);
    int kmax = std::min(m_grid->cell_count_z()-1, idx3D.k + cell_count);

    it.m_index = -1;

    it.m_i_start = imin;
    it.m_j_start = jmin;
    it.m_k_start = kmin;
    it.m_i_end = imax + 1;
    it.m_j_end = jmax + 1;
    it.m_k_end = kmax + 1;

    it.m_i = imin;
    it.m_j = jmin;
    it.m_k = kmin;

    // set iterator to just before the indices of this first cell
    it.m_idx_cell = m_grid->cell_index_1D(Index3D({it.m_i,it.m_j,it.m_k}));
    it.m_it       = m_grid->cell_data()[it.m_idx_cell] - 1;
}

void GridRangeIndexQuery::advance(GridRangeIndexIterator& it)
{
    const auto& cells   = m_grid->cell_data();
    const auto& indices = m_grid->index_data();
    const auto& points  = m_grid->point_data();
    const auto& point   = m_grid->point_data()[m_index];

    while(true)
    {
        // first iterate in current cell
        for(it.m_it = it.m_it + 1; it.m_it < cells[it.m_idx_cell+1]; ++it.m_it)
        {
            int idx_point = indices[it.m_it];

            if(idx_point == m_index) continue;

            if((point - points[idx_point]).squaredNorm() < m_squared_radius)
            {
                it.m_index = idx_point;
                return;
            }
        }

        // advance on next cell
        ++it.m_i;
        if(it.m_i == it.m_i_end)
        {
            it.m_i = it.m_i_start;
            ++it.m_j;
            if(it.m_j == it.m_j_end)
            {
                it.m_j = it.m_j_start;
                ++it.m_k;
                if(it.m_k == it.m_k_end)
                {
                    // end of iteration
                    it.m_index = m_grid->size();
                    return;
                }
            }
        }

        // set iterator to just before the indices of this new cell
        it.m_idx_cell = m_grid->cell_index_1D(Index3D({it.m_i,it.m_j,it.m_k}));
        it.m_it       = m_grid->cell_data()[it.m_idx_cell] - 1;
    }
}

} // namespace pca
