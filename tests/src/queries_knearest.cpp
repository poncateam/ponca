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

/*
template<typename DataPoint>
void testKdTreeKNearestIndex(bool quick = QUICK_TESTS)
{
	using Scalar = typename DataPoint::Scalar;
	using PointContainer = typename KdTreeDense<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	const int k = quick ? 5 : 15;
	auto points = PointContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

    auto kdStart = std::chrono::system_clock::now();
	/// [Kdtree construction and query]
	Ponca::KdTreeDense<DataPoint> kdTree(points);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
        std::vector<int> results; results.reserve( k );
		for (int j : kdTree.k_nearest_neighbors(i, k))
		{
			results.push_back(j);
		}

		bool res = check_k_nearest_neighbors<Scalar>(points, i, k, results);
		VERIFY(res);
	}
    /// [Kdtree construction and query]
    auto kdEnd = std::chrono::system_clock::now();
    auto graphStart = std::chrono::system_clock::now();
    /// [KnnGraph construction and query]
    Ponca::KnnGraph<DataPoint> knnGraph(kdTree, k);
#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        std::vector<int> results; results.reserve( k );
        for (int j : knnGraph.k_nearest_neighbors(i))
        {
            results.push_back(j);
        }

        bool res = check_k_nearest_neighbors<Scalar>(points, i, k, results);
        VERIFY(res);
    }
    /// [KnnGraph construction and query]
    auto graphEnd = std::chrono::system_clock::now();


    std::chrono::duration<double> kdDiff = (kdEnd-kdStart);
    std::chrono::duration<double> graphDiff = (graphEnd-graphStart);

    std::cout << "Test timings: " << "\n"
              << "KdTree   : " <<  kdDiff.count() << "\n"
              << "KnnGraph : " <<  graphDiff.count() << "\n";
}

template<typename DataPoint>
void testKdTreeKNearestPoint(bool quick = QUICK_TESTS)
{
	using Scalar = typename DataPoint::Scalar;
	using PointContainer = typename KdTreeDense<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	const int k = quick ? 5 : 15;
    /// [Kdtree construction]
	auto points = PointContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	KdTreeDense<DataPoint> structure(points);
    /// [Kdtree construction]

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		VectorType point = VectorType::Random();
        std::vector<int> results; results.reserve( k );
		for (int j : structure.k_nearest_neighbors(point, k))
		{
			results.push_back(j);
		}

		bool res = check_k_nearest_neighbors<Scalar>(points, point, k, results);
        VERIFY(res);
	}
}
*/

template<bool doIndexQuery, typename AcceleratingStructure>
void testKdTreeKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const bool quick
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	const int k = quick ? 5 : 15;

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
			return check_k_nearest_neighbors<Scalar>(points, queryInput, k, queryResults);
		}, k
	);
}
template<typename AcceleratingStructure>
void testKnnGraphKNearestNeighbors( AcceleratingStructure& structure,
	typename AcceleratingStructure::PointContainer& points,
	const bool quick
) {
	using DataPoint      = typename AcceleratingStructure::DataPoint;
	using Scalar         = typename DataPoint::Scalar;

	const int k = quick ? 5 : 15;

	testQuery<true, DataPoint>(points,
	[&structure, &k](auto &queryInput) {
			auto mutableQuery = structure.k_nearest_neighbors_empty_index();
			return mutableQuery(queryInput);
		}, [&structure, &k](auto &queryInput) {
			return structure.k_nearest_neighbors(queryInput);
		}, [&points, &k](auto& queryInput, auto& queryResults) {
			return check_k_nearest_neighbors<Scalar>(points, queryInput, k, queryResults);
		}, k
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
	auto kdStart = std::chrono::system_clock::now(); // Only record time for one query
	testKdTreeKNearestNeighbors<true>(kdtreeDense, points, quick);  // Index query test
	auto kdEnd = std::chrono::system_clock::now();
	testKdTreeKNearestNeighbors<false>(kdtreeDense, points, quick); // Position query test

	auto graphStart = std::chrono::system_clock::now();

	const int k = quick ? 5 : 15;
	//////////// Test KnnGraph
	KnnGraph<P> knnGraph(kdtreeDense, k); // We need a large graph, otherwise we might miss some points
											//   (which is the goal of the graph: to replace full euclidean
											//   collection by geodesic-like region growing bounded by
											//   the euclidean ball).

	auto graphEnd = std::chrono::system_clock::now();
	std::chrono::duration<double> kdDiff    = (kdEnd-kdStart);
	std::chrono::duration<double> graphDiff = (graphEnd-graphStart);
	testKnnGraphKNearestNeighbors(knnGraph, points, quick);  // Index query test

	cout << "Timings : " << endl
		 << "    KdTree   : " <<  kdDiff.count() << endl
	     << "    KnnGraph : " <<  graphDiff.count() << endl;
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
