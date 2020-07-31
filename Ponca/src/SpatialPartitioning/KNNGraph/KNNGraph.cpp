#include <PCA/SpacePartitioning/KNNGraph/KNNGraph.h>
#include <PCA/SpacePartitioning/KdTree/KdTree.h>

#include <PCA/Common/Progress.h>

namespace pca {

// KNNGraph --------------------------------------------------------------------

KNNGraph::KNNGraph() :
    m_k(0),
    m_points(nullptr),
    m_indices(nullptr),
    m_verbose(false)
{
}

KNNGraph::KNNGraph(int k) :
    m_k(k),
    m_points(nullptr),
    m_indices(nullptr),
    m_verbose(false)
{
}

void KNNGraph::clear()
{
    m_points  = nullptr;
    m_indices = nullptr;
}

void KNNGraph::build(const KdTree& kdtree)
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

void KNNGraph::build(const KdTree& kdtree, int k)
{
    m_k = k;
    this->build(kdtree);
}

void KNNGraph::build(const KdTree& kdtree, const std::vector<int>& indices)
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

void KNNGraph::build(const KdTree& kdtree, int k, const std::vector<int>& indices)
{
    m_k = k;
    this->build(kdtree, indices);
}

// Query -----------------------------------------------------------------------

KNNGraphQuery KNNGraph::k_nearest_neighbors(int index) const
{
    return KNNGraphQuery(this, index);
}

KNNGraphRangeQuery KNNGraph::range_neighbors(int index, Scalar r) const
{
    return KNNGraphRangeQuery(this, r, index);
}

int KNNGraph::k_neighbor(int idx_point, int i) const
{
    return m_indices->operator[](idx_point * m_k + i);
}

// Empty Query -----------------------------------------------------------------

KNNGraphRangeQuery KNNGraph::range_query(Scalar r) const
{
    return RangeIndexQuery(this, r);
}

// Accessors -------------------------------------------------------------------

int KNNGraph::k() const
{
    return m_k;
}

int KNNGraph::size() const
{
    return m_points->size();
}

const Vector3Array& KNNGraph::point_data() const
{
    return *m_points;
}

Vector3Array& KNNGraph::point_data()
{
    return *m_points;
}

const std::vector<int>& KNNGraph::index_data() const
{
    return *m_indices.get();
}

std::vector<int>& KNNGraph::index_data()
{
    return *m_indices.get();
}

void KNNGraph::set_verbose(bool verbose)
{
    m_verbose = verbose;
}

} // namespace pca
