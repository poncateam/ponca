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

using namespace Ponca;

//! Test kNearestNeighbors query
template<bool doIndexQuery, typename AcceleratingStructure>
auto testKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	std::vector<int>& sample,
	const int retry_number, const int k
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
		}, retry_number, k
	);
}

//! Test kNearestNeighbors query without the k argument (the size of the iterator depends on the acceleration structure (e.g. when using the knnGraph(kdtreeDense, k))
template<typename AcceleratingStructure>
auto testKNearestNeighborsEntirePointSet( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const int retry_number, const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	return testQuery<true, DataPoint>(points,
	[&structure]() {
			return structure.kNearestNeighborsIndexQuery();
		}, [&structure](auto &queryInput) {
			return structure.kNearestNeighbors(queryInput);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return checkKNearestNeighbors<DataPoint>(points, queryInput, k, queryResults);
		}, retry_number
	);
}


template<typename Scalar, int Dim>
void testKNearestNeighborsForAllStructures(const bool quick = QUICK_TESTS)
{
	using P = PointPositionNormal<Scalar, Dim>;
	const int N = quick ? 100 : 1000;
	const int k = quick ? 2 : 15;
	const int retry_number = quick? 1 : 5;
	std::chrono::milliseconds timing;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	cout << endl;
	//////////// Test Dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	timing = testKNearestNeighbors<true>(kdtreeDense, points, sample, retry_number, k);  // Index query test
	cout << "    Compute Time KdTreeDense index query : " <<  timing.count() << "ms" << endl;
	timing = testKNearestNeighbors<false>(kdtreeDense, points, sample, retry_number, k); // Position query test
	cout << "    Compute Time KdTreeDense position query : " <<  timing.count() << "ms" << endl;

	//////////// Test subsample of KdTree
	std::vector<int> subSample;
	KdTreeSparse<P> kdtreeSparse = *buildSubsampledKdTree(points, subSample);
	timing = testKNearestNeighbors<true>(kdtreeSparse, points, subSample, retry_number, k);  // Index query test
	cout << "    Compute Time KdTreeSparse position query : " <<  timing.count() << "ms" << endl;
	timing = testKNearestNeighbors<false>(kdtreeSparse, points, subSample, retry_number, k); // Position query test
	cout << "    Compute Time KdTreeSparse index query : " <<  timing.count() << "ms" << endl;

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, k);
	timing = testKNearestNeighborsEntirePointSet(knnGraph, points, retry_number, k);  // Index query test
	cout << "    Compute Time KnnGraph index query : " <<  timing.count() << "ms" << endl;
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
