/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"
#include "../common/kdtree_utils.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>

using namespace Ponca;

template<typename DataPoint>
void testKdTreeNearestIndex(bool quick = QUICK_TESTS)
{
	using Scalar = typename DataPoint::Scalar;
	using KdTreeType = KdTreeDense<DataPoint>;
	using VectorContainer = typename KdTreeType::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	KdTreeType kdTree(points);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
        std::vector<int> results; results.reserve( 1 );
		for (int j : kdTree.nearest_neighbor(i))
		{
			results.push_back(j);
		}
        VERIFY(results.size() == 1);
		bool res = check_nearest_neighbor<Scalar, VectorContainer>(points, i, results.front());
        VERIFY(res);
	}

    /// [KnnGraph construction]
    Ponca::KnnGraph<DataPoint> knnGraph(kdTree, 1);
    /// [KnnGraph construction]
#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        std::vector<int> results; results.reserve( 1 );
        for (int j : knnGraph.k_nearest_neighbors(i))
        {
            results.push_back(j);
        }
        VERIFY(results.size() == 1);
        bool res = check_nearest_neighbor<Scalar, VectorContainer>(points, i, results.front());
        VERIFY(res);
    }
}

template<typename DataPoint>
void testKdTreeNearestPoint(bool quick = QUICK_TESTS)
{
	using Scalar = typename DataPoint::Scalar;
	using KdTreeType = KdTreeDense<DataPoint>;
	using VectorContainer = typename KdTreeType::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	KdTreeType structure(points);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		VectorType point = VectorType::Random();
        std::vector<int> results; results.reserve( 1 );
		for (int j : structure.nearest_neighbor(point))
		{
			results.push_back(j);
		}
        VERIFY(results.size() == 1);
		bool res = check_nearest_neighbor<Scalar, VectorType, VectorContainer>(points, point, results.front());
        VERIFY(res);
	}
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

    cout << "Test Nearest (from Point) in 3D..." << endl;
	testKdTreeNearestPoint<TestPoint<float, 3>>();
	testKdTreeNearestPoint<TestPoint<double, 3>>();
	testKdTreeNearestPoint<TestPoint<long double, 3>>();

    cout << "Test Nearest (from Point) in 4D..." << endl;
	testKdTreeNearestPoint<TestPoint<float, 4>>();
	testKdTreeNearestPoint<TestPoint<double, 4>>();
	testKdTreeNearestPoint<TestPoint<long double, 4>>();

    cout << "Test Nearest (from Index) in 3D..." << endl;
	testKdTreeNearestIndex<TestPoint<float, 3>>();
	testKdTreeNearestIndex<TestPoint<double, 3>>();
	testKdTreeNearestIndex<TestPoint<long double, 3>>();

    cout << "Test Nearest (from Index) in 4D..." << endl;
	testKdTreeNearestIndex<TestPoint<float, 4>>();
	testKdTreeNearestIndex<TestPoint<double, 4>>();
	testKdTreeNearestIndex<TestPoint<long double, 4>>();
}
