/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../common/testing.h"
#include "../common/kdtree_utils.h"
#include "../split_test_helper.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>
#include <Ponca/src/Common/pointTypes.h>

#define PRINT_TIMING

using namespace Ponca;

template<bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testRangeNeighbors( AcceleratingStructure& structure,
	PointContainer& points,
	std::vector<int>& sample
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	const Scalar r = Eigen::internal::random<Scalar>(Scalar(0.01), Scalar(0.5));

	return testQuery<doIndexQuery, DataPoint>(points,
		[&structure]() {
			if constexpr (doIndexQuery) {
				return structure.rangeNeighborsIndexQuery();
			} else {
				return structure.rangeNeighborsQuery();
			}
		}, [&structure](auto& queryInput, const Scalar _r) {
			return structure.rangeNeighbors(queryInput, _r);
		}, [&points, &sample, &r](auto& queryInput, auto& queryResults) {
			return checkRangeNeighbors<DataPoint>(points, sample, queryInput, r, queryResults);
		}, g_repeat, r
	);
}

template<template <typename> class KdTreeType, typename P>
inline KdTreeType<P> testKdTree(std::vector<P> & points, std::vector<int> & sample) {
	KdTreeType<P> kdtree = *testBuildKdTree<P, KdTreeType>(points, sample);

	std::chrono::milliseconds timing = testRangeNeighbors<true>(kdtree, points, sample);  // Index query test
#ifdef PRINT_TIMING
	cout << "    Compute Time KdTree index query : " <<  timing.count() << "ms" << endl;
#endif
	timing = testRangeNeighbors<false>(kdtree, points, sample); // Position query test
#ifdef PRINT_TIMING
	cout << "    Compute Time KdTree position query : " <<  timing.count() << "ms" << endl;
#endif
	return kdtree;
}

template<typename Scalar, int Dim>
void testRangeNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
	using P = PointPositionNormal<Scalar, Dim>;
	const int N = quick ? 100 : 1000;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	//////////// Test KdTree STL-like containers
	std::vector<int> sampleDense;
	KdTreeDense<P> kdtreeDense = testKdTree<KdTreeDense>(points, sampleDense);
	std::vector<int> sampleSparse;
	testKdTree<KdTreeSparse>(points, sampleSparse);

	//////////// Test KdTree memory array
	std::vector<int> sampleDenseP;
	testKdTree<KdTreeDensePointers>(points, sampleDenseP);
	std::vector<int> sampleSparseP;
	testKdTree<KdTreeSparsePointers>(points, sampleSparseP);

	////////// Test KnnGraph
	 KnnGraph<P> knnGraph(kdtreeDense, N/4); /* We need a large graph, otherwise we might miss some points
	 										   (which is the goal of the graph: to replace full Euclidean
	 										   collection by geodesic-like region growing bounded by
	 										   the Euclidean ball). */
	 std::chrono::milliseconds timing = testRangeNeighbors<true>(knnGraph, points, sampleDense);  // Index query test
 #ifdef PRINT_TIMING
 	cout << "    Compute Time KnnGraph index query : " <<  timing.count() << "ms" << endl;
 #endif
	cout << "  (ok)" << endl;
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
