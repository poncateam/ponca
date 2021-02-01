#pragma once

#include <PCA/SpacePartitioning/Grid/Index3D.h>

#include <PCA/SpacePartitioning/Grid/Query/GridKNearestIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridKNearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridNearestIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridNearestPointQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridRangeIndexQuery.h>
#include <PCA/SpacePartitioning/Grid/Query/GridRangePointQuery.h>

#include <memory>

namespace Ponca {

class Grid
{
    // Types -------------------------------------------------------------------
public:
    using KNearestPointQuery = GridKNearestPointQuery;
    using KNearestIndexQuery = GridKNearestIndexQuery;
    using NearestPointQuery  = GridNearestPointQuery;
    using NearestIndexQuery  = GridNearestIndexQuery;
    using RangePointQuery    = GridRangePointQuery;
    using RangeIndexQuery    = GridRangeIndexQuery;

    // Grid --------------------------------------------------------------------
public:
    Grid();
    Grid(std::shared_ptr<Vector3Array>& points);

    void clear();
    void build(std::shared_ptr<Vector3Array>& points);

    bool valid() const;

    // Query -------------------------------------------------------------------
public:
    KNearestPointQuery k_nearest_neighbors(const Vector3& point, int k) const;
    KNearestIndexQuery k_nearest_neighbors(int index, int k) const;
    NearestPointQuery  nearest_neighbor(const Vector3& point) const;
    NearestIndexQuery  nearest_neighbor(int index) const;
    RangePointQuery    range_neighbors(const Vector3& point, Scalar r) const;
    RangeIndexQuery    range_neighbors(int index, Scalar r) const;

    // Empty Query -------------------------------------------------------------
public:
    KNearestPointQuery k_nearest_point_query(int k = 0) const;
    KNearestIndexQuery k_nearest_index_query(int k = 0) const;
    NearestPointQuery  nearest_point_query() const;
    NearestIndexQuery  nearest_index_query() const;
    RangePointQuery    range_point_query(Scalar r = 0) const;
    RangeIndexQuery    range_index_query(Scalar r = 0) const;

    // Accessors ---------------------------------------------------------------
public:
    int size() const;

    const Vector3Array& point_data() const;
          Vector3Array& point_data();

    const std::vector<int>& cell_data() const;
          std::vector<int>& cell_data();

    const std::vector<int>& index_data() const;
          std::vector<int>& index_data();

    Scalar cell_size() const;
    int cell_count_x() const;
    int cell_count_y() const;
    int cell_count_z() const;

    // Internal ----------------------------------------------------------------
public:
    Index3D cell_index_3D(const Vector3& point) const;
    Index3D cell_index_3D(int idx1D) const;
    int     cell_index_1D(const Vector3& point) const;
    int     cell_index_1D(Index3D idx3D) const;

    Aabb cell_aabb(int idx1D) const;

    // Data --------------------------------------------------------------------
protected:
    std::shared_ptr<Vector3Array>     m_points;
    std::shared_ptr<std::vector<int>> m_cells;
    std::shared_ptr<std::vector<int>> m_indices;

    int m_max_cell_count;

    Aabb   m_aabb;
    Scalar m_cell_size;
    Scalar m_cell_size_inv;
    int    m_cell_count_x;
    int    m_cell_count_y;
    int    m_cell_count_z;
};

} // namespace pca
