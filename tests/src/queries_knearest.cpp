/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/queries_knearest.cpp
 * \brief Test validity of the direct projection on an algebraic sphere
 * \authors Auberval Florian, Nicolas Mellado
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/kdtree_utils.h"
#include "../split_test_helper.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>
#include <Ponca/src/Common/pointTypes.h>

#define PRINT_TIMING

using namespace Ponca;

//! Test kNearestNeighbors query
template <bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testKNearestNeighbors(AcceleratingStructure& structure, PointContainer& points, std::vector<int>& sample,
                           const int k)
{
    using DataPoint = typename AcceleratingStructure::DataPoint;

    return testQuery<doIndexQuery, DataPoint>(
        points,
        [&structure]() {
            if constexpr (doIndexQuery)
            {
                return structure.kNearestNeighborsIndexQuery();
            }
            else
            {
                return structure.kNearestNeighborsQuery();
            }
        },
        [&structure](auto& queryInput, const int _k) { return structure.kNearestNeighbors(queryInput, _k); },
        [&points, &sample, &k](auto& queryInput, auto& queryResults) {
            return checkKNearestNeighbors<DataPoint>(points, sample, queryInput, k, queryResults);
        },
        g_repeat, k);
}

//! \brief Test kNearestNeighbors query without the k argument. The size of the iterator depends on the acceleration
//! structure (e.g. when using the knnGraph(kdtreeDense, k))
template <typename AcceleratingStructure, typename PointContainer>
auto testKNearestNeighborsEntirePointSet(AcceleratingStructure& structure, PointContainer& points, const int k)
{
    using DataPoint = typename AcceleratingStructure::DataPoint;

    return testQuery<true, DataPoint>(
        points, [&structure]() { return structure.kNearestNeighborsIndexQuery(); },
        [&structure](auto& queryInput) { return structure.kNearestNeighbors(queryInput); },
        [&points, &k](auto& queryInput, auto& queryResults) {
            return checkKNearestNeighbors<DataPoint>(points, queryInput, k, queryResults);
        },
        g_repeat);
}

//! \brief Build and test a kdtree for the KNearestNeighbors Query
template <template <typename> class KdTreeType, typename P>
KdTreeType<P> buildAndTestKdTree(std::vector<P>& points, std::vector<int>& sample, const int k,
                                 const std::string& name = "KdTree")
{
    auto kdtree = *testBuildKdTree<P, KdTreeType>(points, sample);

    // Test a kdtree with STL-like vectors
    std::chrono::milliseconds timing = testKNearestNeighbors<true>(kdtree, points, sample, k); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " index query : " << timing.count() << "ms" << endl;
#endif
    timing = testKNearestNeighbors<false>(kdtree, points, sample, k); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " position query : " << timing.count() << "ms" << endl;
#endif

    return kdtree;
}

//! \brief Test a kdtree with pointers for the KNearestNeighbors Query
template <typename P, typename KdTree>
void testStaticKdTree(KdTree& kdtree, const int k, const std::string& name = "KnnGraph")
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

    std::chrono::milliseconds timing =
        testKNearestNeighbors<false>(kdtreeStatic, points, sample, k); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) position query : " << timing.count() << "ms" << endl;
#endif
    timing = testKNearestNeighbors<true>(kdtreeStatic, points, sample, k); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) index query : " << timing.count() << "ms" << endl;
#endif
}

//! \brief Build a KnnGraph and test the KNN query with default container type and raw memory pointers
template <typename P, typename KdTree>
void buildAndTestKnnGraph(KdTree& kdtree, const int k, const std::string& name = "KnnGraph")
{
    auto points = kdtree.points();
    // Test KnnGraph
    KnnGraph<P> knnGraph(kdtree, k);
    std::chrono::milliseconds timing = testKNearestNeighborsEntirePointSet(knnGraph, points, k); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " index query : " << timing.count() << "ms" << endl;
#endif

    // Test the KnnGraph with raw memory pointers
    using KnnGraphPointerStatic = StaticKnnGraphBase<KnnGraphPointerTraits<P>>;
    auto knngraphBuffers        = knnGraph.buffers(); // Buffer that use STL-like containers
    // Convert previous KnnGraph to pointers
    typename KnnGraphPointerStatic::Buffers knnGraphStaticBuffers{
        knngraphBuffers.points.data(), knngraphBuffers.indices.data(), knngraphBuffers.points_size,
        knngraphBuffers.indices_size, k};
    KnnGraphPointerStatic knnGraphStatic(knnGraphStaticBuffers);
    timing = testKNearestNeighborsEntirePointSet(knnGraphStatic, points, k); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time " << name << " (with pointers) index query : " << timing.count() << "ms" << endl;
#endif
}

template <typename Scalar, int Dim>
void testKNearestNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
    using P     = PointPositionNormal<Scalar, Dim>;
    const int N = quick ? 100 : 1000;
    const int k = quick ? 2 : 15;

    //////////// Generate data
    std::vector<P> points(N);
    generateData(points);

    //////////// Test KdTree STL-like containers
    std::vector<int> sampleDense;
    auto kdtreeDense = buildAndTestKdTree<KdTreeDense>(points, sampleDense, k);

    std::vector<int> sampleSparse;
    buildAndTestKdTree<KdTreeSparse>(points, sampleSparse, k, "KdTreeSparse");
    testStaticKdTree<P>(kdtreeDense, k);

    //////////// Test KnnGraph
    buildAndTestKnnGraph<P>(kdtreeDense, k);
    cout << "  (ok)" << endl;
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test kNearestNeighbors query for KdTree and KnnGraph in 3D : " << endl;
    cout << "  float : " << endl;
    CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 3>()));
    cout << "  double : " << endl;
    CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 3>()));
    cout << "  long : " << endl;
    CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 3>()));

    cout << "Test kNearestNeighbors query for KdTree and KnnGraph in 4D : " << endl;
    cout << "  float : " << endl;
    CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 4>()));
    cout << "  double : " << endl;
    CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 4>()));
    cout << "  long : " << endl;
    CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 4>()));

    return EXIT_SUCCESS;
}
