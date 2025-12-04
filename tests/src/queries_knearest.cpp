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

using namespace Ponca;

template<bool doIndexQuery, typename AcceleratingStructure>
void testKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const int retry_number, const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<doIndexQuery, DataPoint>(points,
	[&structure]() {
			if constexpr (doIndexQuery) {
				return structure.kNearestNeighborsIndexQuery();
			} else {
				return structure.kNearestNeighborsQuery();
			}
		}, [&structure](auto &queryInput, const int _k) {
			return structure.kNearestNeighbors(queryInput, _k);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return checkKNearestNeighbors<DataPoint>(points, queryInput, k, queryResults);
		}, retry_number, k
	);
}
template<typename AcceleratingStructure>
void testKNearestNeighborsEntirePointSet( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const int retry_number, const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<true, DataPoint>(points,
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
	using P = TestPoint<Scalar, Dim>;
	const int N = quick ? 100 : 5000;
	const int k = quick ? 2 : 15;
	const int retry_number = quick? 1 : 10;

	//////////// Generate data
	std::vector<P> points(N);
	generateData(points);

	cout << endl;
	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	auto timeStart = std::chrono::system_clock::now(); // Only record time for one query
	testKNearestNeighbors<true>(kdtreeDense, points, retry_number, k);  // Index query test
	cout << "    Compute Time KdTree index query : " <<  (std::chrono::system_clock::now() - timeStart).count() << endl;
	timeStart = std::chrono::system_clock::now();
	testKNearestNeighbors<false>(kdtreeDense, points, retry_number, k); // Position query test
	cout << "    Compute Time KdTree position query : " <<  (std::chrono::system_clock::now() - timeStart).count() << endl;

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, k);
	timeStart = std::chrono::system_clock::now();
	testKNearestNeighborsEntirePointSet(knnGraph, points, retry_number, k);  // Index query test
	cout << "    Compute Time KnnGraph index query : " <<  (std::chrono::system_clock::now() - timeStart).count();
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
		return EXIT_FAILURE;

	cout << "Test kNearestNeighbors query for KdTree and KnnGraph in 3D : " << flush;
	cout << endl << " float :" << flush;
	CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 3>()));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 3>()));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 3>()));

	if (QUICK_TESTS)
		return EXIT_SUCCESS;
	cout << "Test kNearestNeighbors query for KdTree and KnnGraph in 4D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 4>()));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 4>()));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 4>()));

	return EXIT_SUCCESS;
}
