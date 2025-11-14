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

using namespace Ponca;

template<bool doIndexQuery, typename AcceleratingStructure>
void testRangeNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	std::vector<int>& sample
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	Scalar r = Eigen::internal::random<Scalar>(0.01, 0.5);
	
	testQuery<doIndexQuery, DataPoint>(points,
		[&structure, &r](auto &queryInput) {
			if constexpr (doIndexQuery) {
				auto mutableQuery = structure.range_neighbors_empty_index();
				return mutableQuery(queryInput, r);
			} else {
				auto mutableQuery = structure.range_neighbors_empty_position();
				return mutableQuery(queryInput, r);
			}
		}, [&structure, &r](auto &queryInput) {
			return structure.range_neighbors(queryInput, r);
		}, [&points, &sample, &r](auto& queryInput, auto& queryResults) {
			return check_range_neighbors<DataPoint>(points, sample, queryInput, r, queryResults);
		}
	);
}

template<typename Scalar, int Dim>
void testRangeNeighborsForAllStructures(const bool quick)
{
	using P = TestPoint<Scalar, Dim>;
	const int N = quick ? 100 : 5000;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	testRangeNeighbors<true>(kdtreeDense, points, sample);  // Index query test
	testRangeNeighbors<false>(kdtreeDense, points, sample); // Position query test

	//////////// Test subsample of KdTree
	std::vector<int> subSample;
	KdTreeSparse<P> kdtreeSparse = *buildSubsampledKdTree(points, subSample);
	testRangeNeighbors<true>(kdtreeSparse, points, subSample);  // Index query test
	testRangeNeighbors<false>(kdtreeSparse, points, subSample); // Position query test

	// //////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, N/4); /* We need a large graph, otherwise we might miss some points
											   (which is the goal of the graph: to replace full euclidean
											   collection by geodesic-like region growing bounded by
											   the euclidean ball). */
	testRangeNeighbors<true>(knnGraph, points, sample);  // Index query test
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

	cout << "Test range_neighbors query for KdTree and KnnGraph in 3D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testRangeNeighborsForAllStructures<float, 3>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testRangeNeighborsForAllStructures<double, 3>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testRangeNeighborsForAllStructures<long double, 3>(quick)));

	cout << "Test range_neighbors query for KdTree and KnnGraph in 4D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testRangeNeighborsForAllStructures<float, 4>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testRangeNeighborsForAllStructures<double, 4>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testRangeNeighborsForAllStructures<long double, 4>(quick)));
}
