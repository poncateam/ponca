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
auto testRangeNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	std::vector<int>& sample, const int retry_number
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	const Scalar r = Eigen::internal::random<Scalar>(0.01, 0.5);

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
		}, retry_number, r
	);
}

template<typename Scalar, int Dim>
void testRangeNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
	using P = TestPoint<Scalar, Dim>;
	const int N = quick ? 100 : 5000;
	const int retry_number = quick? 1 : 10;
	std::chrono::milliseconds timing;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	auto timeStart = std::chrono::system_clock::now(); // Only record time for one query
	timing = testRangeNeighbors<true>(kdtreeDense, points, sample, retry_number);  // Index query test
	cout << "    Compute Time KdTreeDense index query : " <<  timing.count() << "ms" << endl;
	timing = testRangeNeighbors<false>(kdtreeDense, points, sample, retry_number); // Position query test
	cout << "    Compute Time KdTreeDense position query : " <<  timing.count() << "ms" << endl;

	//////////// Test subsample of KdTree
	std::vector<int> subSample;
	KdTreeSparse<P> kdtreeSparse = *buildSubsampledKdTree(points, subSample);
	timing = testRangeNeighbors<true>(kdtreeSparse, points, subSample, retry_number);  // Index query test
	cout << "    Compute Time KdTreeSparse index query : " <<  timing.count() << "ms" << endl;
	timing = testRangeNeighbors<false>(kdtreeSparse, points, subSample, retry_number); // Position query test
	cout << "    Compute Time KdTreeSparse position query : " <<  timing.count() << "ms" << endl;

	// //////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, N/4); /* We need a large graph, otherwise we might miss some points
											   (which is the goal of the graph: to replace full Euclidean
											   collection by geodesic-like region growing bounded by
											   the Euclidean ball). */
	timing = testRangeNeighbors<true>(knnGraph, points, sample, retry_number);  // Index query test
	cout << "    Compute Time KnnGraph index query : " <<  timing.count() << "ms" << endl;
	cout << "  (ok)" << endl;
}

int main(int argc, char** argv)
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

	if (QUICK_TESTS)
		return EXIT_SUCCESS;

	cout << "Test rangeNeighbors query for KdTree and KnnGraph in 4D : " << endl;
	cout << "  float : " << endl;
	CALL_SUBTEST_1((testRangeNeighborsForAllStructures<float, 4>()));
	cout << "  double : " << endl;
	CALL_SUBTEST_2((testRangeNeighborsForAllStructures<double, 4>()));
	cout << "  long : " << endl;
	CALL_SUBTEST_3((testRangeNeighborsForAllStructures<long double, 4>()));

	return EXIT_SUCCESS;
}
