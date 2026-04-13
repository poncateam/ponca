/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/queries_nearest.cpp
 * \brief Test validity of the direct projection on an algebraic sphere
 * \authors Auberval Florian, Nicolas Mellado
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/kdtree_utils.h"
#include "../split_test_helper.h"

#define PRINT_TIMING

using namespace Ponca;

template <bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testNearestNeighbor(AcceleratingStructure& structure, PointContainer& points, std::vector<int>& sample)
{
    using DataPoint = typename AcceleratingStructure::DataPoint;

    return testQuery<doIndexQuery, DataPoint>(
        points, [&structure](auto& queryInput) { return structure.nearestNeighbor(queryInput); },
        [&points, &sample](auto& queryInput, auto& queryResults) {
            VERIFY((queryResults.size() == 1));
            return checkNearestNeighbor<DataPoint>(points, sample, queryInput, queryResults.front());
        },
        g_repeat);
}

template <bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testFrontOfKNearestNeighbors(AcceleratingStructure& structure, PointContainer& points, std::vector<int>& sample)
{
    using DataPoint = typename AcceleratingStructure::DataPoint;

    return testQuery<doIndexQuery, DataPoint>(
        points, [&structure]() { return structure.kNearestNeighborsIndexQuery(); },
        [&structure](auto& queryInput) { return structure.kNearestNeighbors(queryInput); },
        [&points, &sample](auto& queryInput, auto& queryResults) {
            VERIFY((queryResults.size() == 1));
            return checkNearestNeighbor<DataPoint>(points, sample, queryInput, queryResults.front());
        },
        g_repeat);
}

template <template <typename> class KdTreeType, typename P>
inline KdTreeType<P> testKdTree(std::vector<P>& points, std::vector<int>& sample)
{
    KdTreeType<P> kdtree = *testBuildKdTree<P, KdTreeType>(points, sample);

    std::chrono::milliseconds timing = testNearestNeighbor<true>(kdtree, points, sample); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time KdTree index query : " << timing.count() << "ms" << endl;
#endif
    timing = testNearestNeighbor<false>(kdtree, points, sample); // Position query test
#ifdef PRINT_TIMING
    cout << "    Compute Time KdTree position query : " << timing.count() << "ms" << endl;
#endif
    return kdtree;
}

template <typename Scalar, int Dim>
void testNearestNeighborForAllStructures(const bool quick = QUICK_TESTS)
{
    using P     = PointPositionNormal<Scalar, Dim>;
    const int N = quick ? 100 : 5000;

    //////////// Generate data
    std::vector<P> points(N);
    generateData(points);

    //////////// Test KdTree STL-like containers
    std::vector<int> sampleDense;
    KdTreeDense<P> kdtreeDense = testKdTree<KdTreeDense>(points, sampleDense);
    std::vector<int> sampleSparse;
    testKdTree<KdTreeSparse>(points, sampleSparse);

    //////////// Test KnnGraph
    KnnGraph<P> knnGraph(kdtreeDense, 1);
    std::chrono::milliseconds timing =
        testFrontOfKNearestNeighbors<true>(knnGraph, points, sampleDense); // Index query test
#ifdef PRINT_TIMING
    cout << "    Compute Time KnnGraph index query : " << timing.count() << "ms" << endl;
#endif

    cout << "(ok)";
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test nearestNeighbor query for KdTree and KnnGraph in 3D : " << flush;
    cout << endl << " float : " << flush;
    CALL_SUBTEST_1((testNearestNeighborForAllStructures<float, 3>()));
    cout << endl << " double : " << flush;
    CALL_SUBTEST_2((testNearestNeighborForAllStructures<double, 3>()));
    cout << endl << " long : " << flush;
    CALL_SUBTEST_3((testNearestNeighborForAllStructures<long double, 3>()));

    cout << "Test nearestNeighbor query for KdTree and KnnGraph in 4D : " << flush;
    cout << endl << " float : " << flush;
    CALL_SUBTEST_1((testNearestNeighborForAllStructures<float, 4>()));
    cout << endl << " double : " << flush;
    CALL_SUBTEST_2((testNearestNeighborForAllStructures<double, 4>()));
    cout << endl << " long : " << flush;
    CALL_SUBTEST_3((testNearestNeighborForAllStructures<long double, 4>()));

    return EXIT_SUCCESS;
}
