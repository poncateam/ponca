#include <PCA/SpacePartitioning/Grid/Grid.h>

#include <PCA/Common/Progress.h>
#include <PCA/Common/Assert.h>

#include <list>
#include <numeric>

namespace pca {

// Grid ------------------------------------------------------------------------

Grid::Grid() :
    m_points(nullptr),
    m_cells(nullptr),
    m_indices(nullptr),
    m_max_cell_count(256*256*256),
    m_aabb(),
    m_cell_size(0),
    m_cell_size_inv(0),
    m_cell_count_x(0),
    m_cell_count_y(0),
    m_cell_count_z(0)
{
}

Grid::Grid(std::shared_ptr<Vector3Array>& points) :
    m_points(nullptr),
    m_cells(nullptr),
    m_indices(nullptr),
    m_max_cell_count(256*256*256),
    m_aabb(),
    m_cell_size(0),
    m_cell_size_inv(0),
    m_cell_count_x(0),
    m_cell_count_y(0),
    m_cell_count_z(0)
{
    this->build(points);
}

void Grid::clear()
{
    m_points        = nullptr;
    m_cells         = nullptr;
    m_indices       = nullptr;
    m_aabb.setNull();
    m_cell_size     = 0;
    m_cell_size_inv = 0,
    m_cell_count_x  = 0;
    m_cell_count_y  = 0;
    m_cell_count_z =  0;
}

void Grid::build(std::shared_ptr<Vector3Array>& points)
{
    this->clear();

    m_points = points;

    // compute aabb
    for(const auto& p : *points.get()) m_aabb.extend(p);

    // computes cell size and count
    m_cell_size = std::cbrt(m_aabb.diagonal().prod() / m_max_cell_count);
    m_cell_size_inv = 1./m_cell_size;
    m_cell_count_x = std::ceil(m_aabb.diagonal().x() * m_cell_size_inv);
    m_cell_count_y = std::ceil(m_aabb.diagonal().y() * m_cell_size_inv);
    m_cell_count_z = std::ceil(m_aabb.diagonal().z() * m_cell_size_inv);

    // center aabb
    Vector3 diag   = m_cell_size * Vector3(m_cell_count_x,m_cell_count_y,m_cell_count_z);
    Vector3 center = m_aabb.center();
    m_aabb.min() = center - diag / 2;
    m_aabb.max() = center + diag / 2;

    int cell_count = m_cell_count_x * m_cell_count_y * m_cell_count_z;

    m_indices = std::make_shared<std::vector<int>>(m_points->size());
    std::iota(m_indices->begin(), m_indices->end(), 0);

    m_cells = std::make_shared<std::vector<int>>(cell_count + 1, 0);

    const auto& index_begin = m_indices->begin();
    const auto& cell_begin  = m_cells->begin();
    const auto& cell_end    = m_cells->end();

    auto prog = Progress(m_points->size());
    for(int idx_point=0; idx_point<int(m_points->size()); ++idx_point)
    {
        int idx_cell = this->cell_index_1D(m_points->operator [](idx_point));
        int pos = m_cells->operator [](idx_cell+1);

        std::rotate(index_begin + pos, index_begin + idx_point, index_begin + idx_point + 1);
        std::for_each(cell_begin + idx_cell + 1, cell_end, [](int& i){ ++i ;});

        ++prog;
    }
    PCA_DEBUG_ASSERT(this->valid());
}

bool Grid::valid() const
{
    if(m_points == nullptr)
        return m_cells == nullptr && m_indices == nullptr;

    if(m_cells == nullptr || m_indices == nullptr)
    {
        PCA_DEBUG_ERROR;
        return false;
    }

    if(m_indices->size() != m_points->size())
    {
        PCA_DEBUG_ERROR;
        return false;
    }

    std::vector<bool> b(m_indices->size(), false);
    for(int idx : *m_indices.get())
    {
        if(idx < 0 || int(m_points->size()) <= idx || b[idx])
        {
            PCA_DEBUG_ERROR;
            return false;
        }
        b[idx] = true;
    }
    if(std::find(b.begin(), b.end(), false) != b.end())
    {
        PCA_DEBUG_ERROR;
        return false;
    }

    if(m_cells->back() != int(m_indices->size()))
    {
        PCA_DEBUG_ERROR;
        return false;
    }

    const auto& cells   = *m_cells.get();
    const auto& indices = *m_indices.get();
    const auto& points  = *m_points.get();
    for(int i=0; i<int(cells.size())-1; ++i)
    {
        auto aabb = this->cell_aabb(i);

        for(int j=cells[i]; j<cells[i+1]; ++j)
        {
            int idx_point = indices[j];

            if(!aabb.contains(points[idx_point]))
            {
                PCA_DEBUG_ERROR;
                return false;
            }
        }
    }
    return true;
}

// Query -----------------------------------------------------------------------

GridKNearestPointQuery Grid::k_nearest_neighbors(const Vector3& point, int k) const
{
    return KNearestPointQuery(this, k, point);
}

GridKNearestIndexQuery Grid::k_nearest_neighbors(int index, int k) const
{
    return KNearestIndexQuery(this, k, index);
}

GridNearestPointQuery Grid::nearest_neighbor(const Vector3& point) const
{
    return NearestPointQuery(this, point);
}

GridNearestIndexQuery Grid::nearest_neighbor(int index) const
{
    return NearestIndexQuery(this, index);
}

GridRangePointQuery Grid::range_neighbors(const Vector3& point, Scalar r) const
{
    return RangePointQuery(this, r, point);
}

GridRangeIndexQuery Grid::range_neighbors(int index, Scalar r) const
{
    return RangeIndexQuery(this, r, index);
}

// Empty Query -----------------------------------------------------------------

GridKNearestPointQuery Grid::k_nearest_point_query(int k) const
{
    return KNearestPointQuery(this, k);
}

GridKNearestIndexQuery Grid::k_nearest_index_query(int k) const
{
    return KNearestIndexQuery(this, k);
}

GridNearestPointQuery Grid::nearest_point_query() const
{
    return NearestPointQuery(this);
}

GridNearestIndexQuery Grid::nearest_index_query() const
{
    return NearestIndexQuery(this);
}

GridRangePointQuery Grid::range_point_query(Scalar r) const
{
    return RangePointQuery(this, r);
}

GridRangeIndexQuery Grid::range_index_query(Scalar r) const
{
    return RangeIndexQuery(this, r);
}

// Accessors -------------------------------------------------------------------

int Grid::size() const
{
    return m_points->size();
}

const Vector3Array& Grid::point_data() const
{
    return *m_points.get();
}

Vector3Array& Grid::point_data()
{
    return *m_points.get();
}

const std::vector<int>& Grid::cell_data() const
{
    return *m_cells.get();
}

std::vector<int>& Grid::cell_data()
{
    return *m_cells.get();
}

const std::vector<int>& Grid::index_data() const
{
    return *m_indices.get();
}

std::vector<int>& Grid::index_data()
{
    return *m_indices.get();
}

Scalar Grid::cell_size() const
{
    return m_cell_size;
}

int Grid::cell_count_x() const
{
    return m_cell_count_x;
}

int Grid::cell_count_y() const
{
    return m_cell_count_y;
}

int Grid::cell_count_z() const
{
    return m_cell_count_z;
}

// Internal --------------------------------------------------------------------

Index3D Grid::cell_index_3D(const Vector3& point) const
{
    Index3D idx3D;
    idx3D.i = (point.x() - m_aabb.min().x()) * m_cell_size_inv;
    idx3D.j = (point.y() - m_aabb.min().y()) * m_cell_size_inv;
    idx3D.k = (point.z() - m_aabb.min().z()) * m_cell_size_inv;
    return idx3D;
}

Index3D Grid::cell_index_3D(int idx1D) const
{
    Index3D idx3D;
    idx3D.k = std::floor(idx1D / (m_cell_count_x * m_cell_count_y));
    idx1D  -= idx3D.k * (m_cell_count_x * m_cell_count_y);
    idx3D.j = std::floor(idx1D / m_cell_count_x);
    idx1D  -= idx3D.j * m_cell_count_x;
    idx3D.i = idx1D;
    return idx3D;
}

int Grid::cell_index_1D(const Vector3& point) const
{
    return this->cell_index_1D(this->cell_index_3D(point));
}

int Grid::cell_index_1D(Index3D idx3D) const
{
    return idx3D.k * (m_cell_count_x * m_cell_count_y) + idx3D.j * m_cell_count_x + idx3D.i;
}

Aabb Grid::cell_aabb(int idx1D) const
{
    auto idx3D = this->cell_index_3D(idx1D);

    Aabb aabb;
    aabb.min().x() = m_aabb.min().x() + m_cell_size * idx3D.i;
    aabb.min().y() = m_aabb.min().y() + m_cell_size * idx3D.j;
    aabb.min().z() = m_aabb.min().z() + m_cell_size * idx3D.k;

    aabb.max().x() = m_aabb.min().x() + m_cell_size * (idx3D.i + 1);
    aabb.max().y() = m_aabb.min().y() + m_cell_size * (idx3D.j + 1);
    aabb.max().z() = m_aabb.min().z() + m_cell_size * (idx3D.k + 1);

    return aabb;
}

} // namespace pca
