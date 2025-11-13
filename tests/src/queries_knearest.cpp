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
	const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<doIndexQuery, DataPoint>(points,
	[&structure, &k](auto &queryInput) {
			if constexpr (doIndexQuery) {
				auto mutableQuery = structure.k_nearest_neighbors_empty_index();
				return mutableQuery(queryInput, k);
			} else {
				auto mutableQuery = structure.k_nearest_neighbors_empty_position();
				return mutableQuery(queryInput, k);
			}
		}, [&structure, &k](auto &queryInput) {
			return structure.k_nearest_neighbors(queryInput, k);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return check_k_nearest_neighbors<DataPoint>(points, queryInput, k, queryResults);
		}, k
	);
}
template<typename AcceleratingStructure>
void testKNearestNeighborsEntirePointSet( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const int k
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;

	testQuery<true, DataPoint>(points,
	[&structure](auto &queryInput) {
			auto mutableQuery = structure.k_nearest_neighbors_empty_index();
			return mutableQuery(queryInput);
		}, [&structure](auto &queryInput) {
			return structure.k_nearest_neighbors(queryInput);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return check_k_nearest_neighbors<DataPoint>(points, queryInput, k, queryResults);
		}, k
	);
}


template<typename Scalar, int Dim>
void testKNearestNeighborsForAllStructures(const bool quick)
{
	using P = TestPoint<Scalar, Dim>;
	// Generate data
	const int N = quick ? 100 : 5000;
	const int k = quick ? 5 : 15;

	std::vector<P> points(N);
	generateData(points);

	cout << endl;
	//////////// Test dense KdTree
	std::vector<int> sample;
	KdTreeDense<P> kdtreeDense = *buildKdTreeDense<P>(points, sample);
	auto timeStart = std::chrono::system_clock::now(); // Only record time for one query
	testKNearestNeighbors<true>(kdtreeDense, points, k);  // Index query test
	cout << "    Compute Time KdTree index query : " <<  (std::chrono::system_clock::now() - timeStart).count() << endl;
	timeStart = std::chrono::system_clock::now();
	testKNearestNeighbors<false>(kdtreeDense, points, k); // Position query test
	cout << "    Compute Time KdTree position query : " <<  (std::chrono::system_clock::now() - timeStart).count() << endl;

	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, k);
	timeStart = std::chrono::system_clock::now();
	testKNearestNeighborsEntirePointSet(knnGraph, points, k);  // Index query test
	cout << "    Compute Time KnnGraph index query : " <<  (std::chrono::system_clock::now() - timeStart).count();
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
		return EXIT_FAILURE;

#ifndef NDEBUG
	bool quick = true;
#else
	bool quick = false;
#endif

	cout << "Test k_nearest_neighbors query for KdTree and KnnGraph in 3D : " << flush;
	cout << endl << " float :" << flush;
	CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 3>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 3>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 3>(quick)));

	cout << "Test k_nearest_neighbors query for KdTree and KnnGraph in 4D : " << flush;
	cout << endl << " float : " << flush;
	CALL_SUBTEST_1((testKNearestNeighborsForAllStructures<float, 4>(quick)));
	cout << endl << " double : " << flush;
	CALL_SUBTEST_2((testKNearestNeighborsForAllStructures<double, 4>(quick)));
	cout << endl << " long : " << flush;
	CALL_SUBTEST_3((testKNearestNeighborsForAllStructures<long double, 4>(quick)));
}
