/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "knnGraph.h"
#include "../KdTree/KdTree.h"

namespace Ponca {

// knnGraph --------------------------------------------------------------------

KnnGraph::KnnGraph() :
    m_k(0),
    m_points(nullptr),
    m_indices(nullptr),
    m_verbose(false)
{
}

KnnGraph::KnnGraph(int k) :
    m_k(k),
    m_points(nullptr),
    m_indices(nullptr),
    m_verbose(false)
{
}

void KnnGraph::clear()
{
    m_points  = nullptr;
    m_indices = nullptr;
}

void KnnGraph::build(const KdTree& kdtree)
{
    this->clear();

    const int size = kdtree.size();

    m_points = kdtree.point_ptr();

    m_indices = std::make_shared<std::vector<int>>(size * m_k, -1);
    auto& indices = *m_indices.get();

    auto q = kdtree.k_nearest_index_query(m_k);

    auto prog = Progress(size, m_verbose);
    #pragma omp parallel for firstprivate(q)
    for(int i=0; i<size; ++i)
    {
        q.set_index(i);

        int j = 0;
        for(int n : q)
        {
            indices[i * m_k + j] = n;
            ++j;
        }
        ++prog;
    }
}

void KnnGraph::build(const KdTree& kdtree, int k)
{
    m_k = k;
    this->build(kdtree);
}

void KnnGraph::build(const KdTree& kdtree, const std::vector<int>& indices)
{
    this->clear();

    const int size = indices.size();

    m_points = kdtree.point_ptr();

    m_indices = std::make_shared<std::vector<int>>(size * m_k, -1);

    auto q = kdtree.k_nearest_index_query(m_k);

    auto prog = Progress(size, m_verbose);
    #pragma omp parallel for firstprivate(q)
    for(int i=0; i<size; ++i)
    {
        q.set_index(indices[i]);

        int j = 0;
        for(int n : q)
        {
            (*m_indices)[i * m_k + j] = n;
            ++j;
        }
        ++prog;
    }
}

void KnnGraph::build(const KdTree& kdtree, int k, const std::vector<int>& indices)
{
    m_k = k;
    this->build(kdtree, indices);
}

// Query -----------------------------------------------------------------------

knnGraphQuery KnnGraph::k_nearest_neighbors(int index) const
{
    return KnnGraphQuery(this, index);
}

knnGraphRangeQuery KnnGraph::range_neighbors(int index, Scalar r) const
{
    return knnGraphRangeQuery(this, r, index);
}

int KnnGraph::k_neighbor(int idx_point, int i) const
{
    return m_indices->operator[](idx_point * m_k + i);
}

// Empty Query -----------------------------------------------------------------

knnGraphRangeQuery KnnGraph::range_query(Scalar r) const
{
    return RangeIndexQuery(this, r);
}

// Accessors -------------------------------------------------------------------

int KnnGraph::k() const
{
    return m_k;
}

int KnnGraph::size() const
{
    return m_points->size();
}

const Vector3Array& KnnGraph::point_data() const
{
    return *m_points;
}

Vector3Array& KnnGraph::point_data()
{
    return *m_points;
}

const std::vector<int>& KnnGraph::index_data() const
{
    return *m_indices.get();
}

std::vector<int>& KnnGraph::index_data()
{
    return *m_indices.get();
}

void KnnGraph::set_verbose(bool verbose)
{
    m_verbose = verbose;
}

} // namespace Ponca
