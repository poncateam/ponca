#pragma once

#include "Query/knnGraphQuery.h"
#include "Query/knnGraphRangeQuery.h"

#include <memory>

namespace Ponca {

class KdTree;

class KnnGraph
{
    // Types -------------------------------------------------------------------
public:
    using KNearestIndexQuery = KnnGraphQuery;
    using RangeIndexQuery    = KnnGraphRangeQuery;

    // knnGraph ----------------------------------------------------------------
public:
    KnnGraph();
    KnnGraph(int k);

    void clear();
    void build(const KdTree& kdtree);
    void build(const KdTree& kdtree, int k);
    void build(const KdTree& kdtree, const std::vector<int>& indices);
    void build(const KdTree& kdtree, int k, const std::vector<int>& indices);


    // Query -------------------------------------------------------------------
public:
    KNearestIndexQuery k_nearest_neighbors(int index) const;
    RangeIndexQuery    range_neighbors(int index, Scalar r) const;

    int k_neighbor(int index, int i) const;

    // Empty Query -------------------------------------------------------------
public:
    RangeIndexQuery range_query(Scalar r = 0) const;

    // Accessors ---------------------------------------------------------------
public:
    int k() const;
    int size() const;

    const Vector3Array& point_data() const;
          Vector3Array& point_data();

    const std::vector<int>& index_data() const;
          std::vector<int>& index_data();

    void set_verbose(bool verbose = true);

    // Data --------------------------------------------------------------------
protected:
    int m_k;
    std::shared_ptr<Vector3Array>     m_points;
    std::shared_ptr<std::vector<int>> m_indices;

    bool m_verbose;
};

} // namespace Ponca
