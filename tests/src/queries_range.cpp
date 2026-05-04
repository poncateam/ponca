/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/queries_range.cpp
 * \brief Test validity of the direct projection on an algebraic sphere
 * \authors Auberval Florian, Nicolas Mellado
 */

#include "../common/testing.h"
#include "../common/kdtree_utils.h"
#include "../split_test_helper.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>
#include <Ponca/src/Common/pointTypes.h>

#define PRINT_TIMING

using namespace Ponca;

template <bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testRangeNeighbors(AcceleratingStructure& structure, PointContainer& points, std::vector<int>& sample)
{
    using DataPoint = typename AcceleratingStructure::DataPoint;
    using Scalar    = typename DataPoint::Scalar;

    const Scalar r = Eigen::internal::random<Scalar>(Scalar(0.01), Scalar(0.5));

    return testQuery<doIndexQuery, DataPoint>(
        points,
        [&structure]() {
            if constexpr (doIndexQuery)
            {
                return structure.rangeNeighborsIndexQuery();
            }
            else
            {
                return structure.rangeNeighborsQuery();
            }
        },
        [&structure](auto& queryInput, const Scalar _r) { return structure.rangeNeighbors(queryInput, _r); },
        [&points, &sample, &r](auto& queryInput, auto& queryResults) {
            return checkRangeNeighbors<DataPoint>(points, sample, queryInput, r, queryResults);
        },
        g_repeat, r);
}

//! \brief Build and test a kdtree for the rangeNeighbors Query
template <template <typename> class KdTreeType, typename P>
KdTreeType<P> buildAndTestKdTree(std::vector<P>& points, std::vector<int>& sample, const std::string& name = "KdTree")
{
    KdTreeType<P> kdtree = std::move(*testBuildKdTree<P, KdTreeType>(points, sample));

    std::chrono::milliseconds timing = testRangeNeighbors<true>(kdtree, points, sample); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " index query : " << timing.count() << "ms" << endl;
#endif
    timing = testRangeNeighbors<false>(kdtree, points, sample); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " position query : " << timing.count() << "ms" << endl;
#endif
    return kdtree;
}

//! \brief Test a kdtree with pointers for the KNearestNeighbors Query
template <typename P, typename KdTree>
void testStaticKdTree(KdTree& kdtree, const std::string& name = "KdTree")
{
    auto points = kdtree.points();
    auto sample = kdtree.samples();

    // Test the KdTree with raw memory pointers
    using KdTreePointerStatic = StaticKdTreeBase<KdTreePointerTraits<P>>;
    auto kdtreeBuffers        = kdtree.buffers(); // Buffer that use STL-like containers
    // Convert previous KdTree to pointers
    typename KdTreePointerStatic::Buffers kdtreeStaticBuffers{kdtreeBuffers.points.data(),  kdtreeBuffers.nodes.data(),
                                                              kdtreeBuffers.indices.data(), kdtreeBuffers.points_size,
                                                              kdtreeBuffers.nodes_size,     kdtreeBuffers.indices_size};
    KdTreePointerStatic kdtreeStatic(kdtreeStaticBuffers);

    std::chrono::milliseconds timing = testRangeNeighbors<false>(kdtreeStatic, points, sample); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) position query : " << timing.count() << "ms" << endl;
#endif
    timing = testRangeNeighbors<true>(kdtreeStatic, points, sample); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) index query : " << timing.count() << "ms" << endl;
#endif
}

//! \brief Build a KnnGraph and test the KNN query with default container type and raw memory pointers
template <typename P, typename KdTree>
void buildAndTestKnnGraph(KdTree& kdtree, std::vector<int>& sampleDense, const std::string& name = "KnnGraph")
{
    auto points = kdtree.points();
    const int k = std::min(100, kdtree.pointCount() / 4);
    // Test KnnGraph
    KnnGraph<P> knnGraph(kdtree, k); /* We need a large graph, otherwise we might miss some points
                                        (which is the goal of the graph: to replace full Euclidean
                                        collection by geodesic-like region growing bounded by
                                        the Euclidean ball). */
    std::chrono::milliseconds timing = testRangeNeighbors<true>(knnGraph, points, sampleDense); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " index query : " << timing.count() << "ms" << endl;
#endif
    cout << "  (ok)" << endl;

    // Test the KnnGraph with raw memory pointers
    using KnnGraphPointerStatic = StaticKnnGraphBase<KnnGraphPointerTraits<P>>;
    auto knngraphBuffers        = knnGraph.buffers(); // Buffer that use STL-like containers
    // Convert previous KnnGraph to pointers
    const P* pts = knngraphBuffers.points.data();
    typename KnnGraphPointerStatic::Buffers knnGraphStaticBuffers{
        pts, knngraphBuffers.indices.data(), knngraphBuffers.points_size,
        knngraphBuffers.indices_size, k};
    KnnGraphPointerStatic knnGraphStatic(knnGraphStaticBuffers);
    timing = testRangeNeighbors<true>(knnGraphStatic, points, sampleDense); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) index query : " << timing.count() << "ms" << endl;
#endif
    cout << "  (ok)" << endl;
}

template <typename Scalar, int Dim>
void testRangeNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
    using P     = PointPositionNormal<Scalar, Dim>;
    const int N = quick ? 100 : 1000;

    //////////// Generate data
    std::vector<P> points(N);
    generateData(points);

    //////////// Test KdTree STL-like containers
    std::vector<int> sampleDense;
    auto kdtreeDense = buildAndTestKdTree<KdTreeDense>(points, sampleDense);

    std::vector<int> sampleSparse;
    buildAndTestKdTree<KdTreeSparse>(points, sampleSparse, "KdTreeSparse");
    testStaticKdTree<P>(kdtreeDense);

    ////////// Test KnnGraph
    buildAndTestKnnGraph<P>(kdtreeDense, sampleDense);
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test rangeNeighbors query for KdTree and KnnGraph in 3D : " << endl;
    cout << "  float : " << endl;
    CALL_SUBTEST_1((testRangeNeighborsForAllStructures<float, 3>()));
    cout << "  double : " << endl;
    CALL_SUBTEST_2((testRangeNeighborsForAllStructures<double, 3>()));
    cout << "  long : " << endl;
    CALL_SUBTEST_3((testRangeNeighborsForAllStructures<long double, 3>()));

    cout << "Test rangeNeighbors query for KdTree and KnnGraph in 4D : " << endl;
    cout << "  float : " << endl;
    CALL_SUBTEST_1((testRangeNeighborsForAllStructures<float, 4>()));
    cout << "  double : " << endl;
    CALL_SUBTEST_2((testRangeNeighborsForAllStructures<double, 4>()));
    cout << "  long : " << endl;
    CALL_SUBTEST_3((testRangeNeighborsForAllStructures<long double, 4>()));

    return EXIT_SUCCESS;
}
