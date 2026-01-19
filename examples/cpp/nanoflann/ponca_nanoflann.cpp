

/*!
\file examples/cpp/nanoflann/ponca_nanoflann.cpp
\brief Demonstrate how to use Ponca with Nanoflann, and compare performances with built-in kdtree


\author: Nicolas Mellado
*/

#include <iostream>
#include <chrono>
#include <vector>

#include "./nanoflann.hpp"

#include <Ponca/Fitting>
#include <Ponca/SpatialPartitioning>


using namespace Ponca;


// This class defines the input data format
class MyPoint
{
public:
    enum {Dim = 3};
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

    PONCA_MULTIARCH inline MyPoint(const VectorType& _pos    = VectorType::Zero()) : m_pos(_pos) {}
    PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }

    // Generate points on a sphere of radius 1
    static inline MyPoint Random()
    {
        return { VectorType::Random().normalized() * Eigen::internal::random<Scalar>(0.9,1.1) };
    }

private:
    VectorType m_pos;
};

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

//! [Define Fit Type]
using NeighborFilter =  DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<Scalar> > ;
using FitType    = Basket<MyPoint,NeighborFilter, Ponca::DryFit>; // build a fitting object that does nothing
//! [Define Fit Type]


///// Nanoflann stuff
struct NFPointCloud
{
    using coord_t = typename MyPoint::Scalar;  //!< The type of each coordinate

    const std::vector<MyPoint>& pts;

    inline NFPointCloud(const std::vector<MyPoint>& _pts) : pts(_pts) {}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return pts[idx].pos()[0];
        else if (dim == 1)
            return pts[idx].pos()[1];
        else
            return pts[idx].pos()[2];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor< Scalar, NFPointCloud>, NFPointCloud, 3>;
///// Enf of nanoflann stuff

int test_raw(FitType& f, const std::vector<MyPoint>& _vecs)
{
    if(! (f.compute( _vecs ) == STABLE) )
        std::cerr << "[raw] Something weird happened" << std::endl;
    return f.getNumNeighbors();
}

int test_ponca_kdtree(FitType& f, const std::vector<MyPoint>& _vecs, VectorType _p, const KdTree<MyPoint>& tree, Scalar tmax){
    if(! (
            //! [Use Ponca KdTree]
f.computeWithIds( tree.rangeNeighbors(_p, tmax), _vecs )
            //! [Use Ponca KdTree]
            == STABLE) )
        std::cerr << "[ponca_kdtree] Something weird happened" << std::endl;
    return f.getNumNeighbors();
}

int test_nanflann_kdtree(FitType& f, const std::vector<MyPoint>& _vecs, VectorType _p, const my_kd_tree_t& tree, Scalar tmax)
{
    // radius search:
    const Scalar                                       squaredRadius = 1;
    std::vector<nanoflann::ResultItem<size_t, Scalar>> indices_dists;
    nanoflann::RadiusResultSet<Scalar, size_t>         resultSet(
            tmax*tmax, indices_dists);

    tree.findNeighbors(resultSet, _p.data());

    //! [Use NanoFlann KdTree]
f.init();
// Compute
auto res = Ponca::UNDEFINED;
do {
    f.startNewPass();
    for (const auto& r : resultSet.m_indices_dists){
        f.addNeighbor(_vecs[r.first]);
    }
    res = f.finalize();
} while ( res == NEED_OTHER_PASS );
    //! [Use NanoFlann KdTree]

    if(res != STABLE)
        std::cerr << "[nanoflann_kdtree] Something weird happened" << std::endl;
    return f.getNumNeighbors();
}


int main()
{
    using Scalar = typename MyPoint::Scalar;
    // init random point cloud
    int n = 100000;
    std::vector<MyPoint> vecs (n);
    std::generate(vecs.begin(), vecs.end(), []() {return MyPoint::Random(); });

    //! [Create Ponca KdTree]
KdTreeDense<MyPoint> ponca_tree(vecs);
    //! [Create Ponca KdTree]

    //! [Create NanoFlann KdTree]
NFPointCloud nfcloud(vecs);
my_kd_tree_t mat_index(3, nfcloud);
    //! [Create NanoFlann KdTree]

    Scalar tmax = 0.2;

    FitType fit;

    int nbrun = 1000;
    std::vector<typename MyPoint::VectorType> queries (nbrun);
    std::generate(queries.begin(), queries.end(), []() {return MyPoint::Random().pos(); });

    int neiRaw {0}, neiPonca {0}, neiFlann {0};
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i != nbrun; ++i)
    {
        fit.setNeighborFilter({queries[i], tmax});
        neiRaw += test_raw(fit, vecs);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> rawDiff = (end-start);

    start = std::chrono::system_clock::now();
    for(int i = 0; i != nbrun; ++i)
    {
        fit.setNeighborFilter({queries[i], tmax});
        neiPonca += test_ponca_kdtree(fit, vecs, queries[i], ponca_tree, tmax);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> poncaDiff = (end-start);

    start = std::chrono::system_clock::now();
    for(int i = 0; i != nbrun; ++i)
    {
        fit.setNeighborFilter({queries[i], tmax});
        neiFlann += test_nanflann_kdtree(fit, vecs, queries[i], mat_index, tmax);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> nanoflannDiff = (end-start);

    std::cout << "Timings: " << "\n"
              << "Raw :       " <<  rawDiff.count() << "\n"
              << "Ponca :     " <<  poncaDiff.count() << "\n"
              << "Nanoflann : " <<  nanoflannDiff.count() << "\n";
    std::cout << "Number of neighbors: " << "\n"
              << "Raw :       " <<  neiRaw << "\n"
              << "Ponca :     " <<  neiPonca << "\n"
              << "Nanoflann : " <<  neiFlann << "\n";
}
