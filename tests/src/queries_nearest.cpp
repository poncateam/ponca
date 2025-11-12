/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"
#include "../common/kdtree_utils.h"
#include "../split_test_helper.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>

using namespace Ponca;

template<bool doIndexQuery, typename AcceleratingStructure>
void testNearestNeighbor( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	testQuery<doIndexQuery, DataPoint>(points, [&structure](auto &queryInput) {
			return structure.nearest_neighbor(queryInput);
		}, [&points](auto& queryInput, auto& queryResults) {
			VERIFY((queryResults.size() == 1));
			return check_nearest_neighbor<Scalar>(points, queryInput, queryResults.front());
		}, 1
	);
}

template<bool doIndexQuery, typename AcceleratingStructure>
void testFrontOfKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	testQuery<doIndexQuery, DataPoint>(points, [&structure](auto &queryInput) {
			return structure.k_nearest_neighbors(queryInput);
		}, [&points](auto& queryInput, auto& queryResults) {
			return check_nearest_neighbor<Scalar>(points, queryInput, queryResults.front());
		}, 1
	);
}

template<typename Scalar, int Dim>
void testNearestNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
	using P = TestPoint<Scalar, Dim>;

	// Generate data
	const int N = quick ? 100 : 500;
	std::vector<P> points(N);
	generateData(points);

	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	testNearestNeighbor<true>(kdtreeDense, points);  // Index query test
	testNearestNeighbor<false>(kdtreeDense, points); // Position query test

	//////////// Test subsample of KdTree
	// std::vector<int> subSample;
	// KdTreeSparse<P> kdtreeSparse = *buildSubsampledKdTree(points, subSample);
	// testNearestNeighbors<true>(kdtreeSparse, points, subSample);  // Index query test
	// testNearestNeighbors<false>(kdtreeSparse, points, subSample); // Position query test

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, 1);
	testFrontOfKNearestNeighbors<true>(knnGraph, points);  // Index query test
}
int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

	cout << "Test range_neighbors query for KdTree and KnnGraph in 3D : " << endl;
	cout << " float" << flush;
	CALL_SUBTEST_1((testNearestNeighborsForAllStructures<float, 3>()));
	cout << " (ok), double" << flush;
	CALL_SUBTEST_2((testNearestNeighborsForAllStructures<double, 3>()));
	cout << " (ok), long " << flush;
	CALL_SUBTEST_3((testNearestNeighborsForAllStructures<long double, 3>()));
	cout << " (ok)." << flush << endl;

	if (QUICK_TESTS)
		return EXIT_SUCCESS;

	cout << "Test range_neighbors query for KdTree and KnnGraph in 4D : " << endl;
	cout << " float" << flush;
	CALL_SUBTEST_1((testNearestNeighborsForAllStructures<float, 4>()));
	cout << " (ok), double" << flush;
	CALL_SUBTEST_2((testNearestNeighborsForAllStructures<double, 4>()));
	cout << " (ok), long " << flush;
	CALL_SUBTEST_3((testNearestNeighborsForAllStructures<long double, 4>()));
	cout << " (ok)." << flush << endl;

	return EXIT_SUCCESS;
}
