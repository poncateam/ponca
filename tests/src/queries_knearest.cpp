/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
template<bool doIndexQuery, typename AcceleratingStructure, typename PointContainer>
auto testKNearestNeighbors( AcceleratingStructure& structure,
	PointContainer& points,
	std::vector<int>& sample,
	const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	return testQuery<doIndexQuery, DataPoint>(points,
	[&structure]() {
			if constexpr (doIndexQuery) {
				return structure.kNearestNeighborsIndexQuery();
			} else {
				return structure.kNearestNeighborsQuery();
			}
		}, [&structure](auto &queryInput, const int _k) {
			return structure.kNearestNeighbors(queryInput, _k);
		}, [&points, &sample, &k](auto& queryInput, auto& queryResults) {
			return checkKNearestNeighbors<DataPoint>(points, sample, queryInput, k, queryResults);
		}, g_repeat, k
	);
}

//! Test kNearestNeighbors query without the k argument (the size of the iterator depends on the acceleration structure (e.g. when using the knnGraph(kdtreeDense, k))
template<typename AcceleratingStructure, typename PointContainer>
auto testKNearestNeighborsEntirePointSet( AcceleratingStructure& structure,
	PointContainer& points,
	const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	return testQuery<true, DataPoint>(points,
	[&structure]() {
			return structure.kNearestNeighborsIndexQuery();
		}, [&structure](auto &queryInput) {
			return structure.kNearestNeighbors(queryInput);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return checkKNearestNeighbors<DataPoint>(points, queryInput, k, queryResults);
		}, g_repeat
	);
}

template<template <typename> class KdTreeType, typename P>
inline KdTreeType<P> testKdTree(std::vector<P> & points, std::vector<int> & sample, const int k) {
	KdTreeType<P> kdtree = *testBuildKdTree<P, KdTreeType>(points, sample);

	std::chrono::milliseconds timing = testKNearestNeighbors<true>(kdtree, points, sample, k);  // Index query test
#ifdef PRINT_TIMING
	cout << "    Compute Time KdTree index query : " <<  timing.count() << "ms" << endl;
#endif
	timing = testKNearestNeighbors<false>(kdtree, points, sample, k); // Position query test
#ifdef PRINT_TIMING
	cout << "    Compute Time KdTree position query : " <<  timing.count() << "ms" << endl;
#endif
	return kdtree;
}

template<typename Scalar, int Dim>
void testKNearestNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
	using P = PointPositionNormal<Scalar, Dim>;
	const int N = quick ? 100 : 1000;
	const int k = quick ? 2 : 15;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	cout << endl;
	//////////// Test KdTree STL-like containers
	std::vector<int> sampleDense;
	KdTreeDense<P> kdtreeDense = testKdTree<KdTreeDense>(points, sampleDense, k);
	std::vector<int> sampleSparse;
	testKdTree<KdTreeSparse>(points, sampleSparse, k);

	//////////// Test KdTree memory array
	std::vector<int> sampleDenseP;
	testKdTree<KdTreeDensePointers>(points, sampleDenseP, k);
	std::vector<int> sampleSparseP;
	testKdTree<KdTreeSparsePointers>(points, sampleSparseP, k);

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, k);

	std::chrono::milliseconds timing = testKNearestNeighborsEntirePointSet(knnGraph, points, k);  // Index query test
#ifdef PRINT_TIMING
	cout << "    Compute Time KnnGraph index query : " <<  timing.count() << "ms" << endl;
#endif
	cout << "  (ok)" << endl;
}

int main(int argc, char** argv)
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
