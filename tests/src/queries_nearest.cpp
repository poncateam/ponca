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
	typename AcceleratingStructure::PointContainer& points,
	std::vector<int>& sample
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<doIndexQuery, DataPoint>(points, [&structure](auto &queryInput) {
			return structure.nearest_neighbor(queryInput);
		}, [&points, &sample](auto& queryInput, auto& queryResults) {
			VERIFY((queryResults.size() == 1));
			return check_nearest_neighbor<DataPoint>(points, sample, queryInput, queryResults.front());
		}, 1
	);
}

template<bool doIndexQuery, typename AcceleratingStructure>
void testFrontOfKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	std::vector<int>& sample
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<doIndexQuery, DataPoint>(points, [&structure](auto &queryInput) {
			return structure.k_nearest_neighbors(queryInput);
		}, [&points, &sample](auto& queryInput, auto& queryResults) {
			return check_nearest_neighbor<DataPoint>(points, sample, queryInput, queryResults.front());
		}, 1
	);
}

template<typename Scalar, int Dim>
void testNearestNeighborForAllStructures(const bool quick)
{
	using P = TestPoint<Scalar, Dim>;
	const int N = quick ? 100 : 5000;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	testNearestNeighbor<true>(kdtreeDense, points, sample);  // Index query test
	testNearestNeighbor<false>(kdtreeDense, points, sample); // Position query test

	//////////// Test subsample of KdTree
	std::vector<int> subSample;
	KdTreeSparse<P> kdtreeSparse = *buildSubsampledKdTree(points, subSample);
	testNearestNeighbor<true>(kdtreeSparse, points, subSample);  // Index query test
	testNearestNeighbor<false>(kdtreeSparse, points, subSample); // Position query test

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, 1);
	testFrontOfKNearestNeighbors<true>(knnGraph, points, sample);  // Index query test

	cout << "(ok)";
}
int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

#ifndef NDEBUG
    const bool quick = true;
#else
    const bool quick = false;
#endif

	cout << "Test nearest_neighbor query for KdTree and KnnGraph in 3D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testNearestNeighborForAllStructures<float, 3>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testNearestNeighborForAllStructures<double, 3>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testNearestNeighborForAllStructures<long double, 3>(quick)));

	cout << "Test nearest_neighbor query for KdTree and KnnGraph in 4D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testNearestNeighborForAllStructures<float, 4>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testNearestNeighborForAllStructures<double, 4>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testNearestNeighborForAllStructures<long double, 4>(quick)));
}
