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

using namespace Ponca;

template<typename DataPoint>
void testKdTreeKNearestIndex(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	const int k = quick ? 5 : 15;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	/// [Kdtree construction and query]
	Ponca::KdTree<DataPoint> structure(points);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
        std::vector<int> results; results.reserve( k );
		for (int j : structure.k_nearest_neighbors(i, k))
		{
			results.push_back(j);
		}

		bool res = check_k_nearest_neighbors<Scalar, VectorContainer>(points, i, k, results);
		VERIFY(res);
	}
    /// [Kdtree construction and query]
}

template<typename DataPoint>
void testKdTreeKNearestPoint(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	const int k = quick ? 5 : 15;
    /// [Kdtree construction]
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	KdTree<DataPoint> structure(points);
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

		bool res = check_k_nearest_neighbors<Scalar, VectorType, VectorContainer>(points, point, k, results);
        VERIFY(res);
	}
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

    cout << "Test KNearest (from Point) in 3D..." << endl;
	testKdTreeKNearestPoint<TestPoint<float, 3>>(false);
	testKdTreeKNearestPoint<TestPoint<double, 3>>(false);
	testKdTreeKNearestPoint<TestPoint<long double, 3>>(false);

    cout << "Test KNearest (from Point) in 3D..." << endl;
	testKdTreeKNearestPoint<TestPoint<float, 4>>(false);
	testKdTreeKNearestPoint<TestPoint<double, 4>>(false);
	testKdTreeKNearestPoint<TestPoint<long double, 4>>(false);

    cout << "Test KNearest (from Index) in 3D..." << endl;
	testKdTreeKNearestIndex<TestPoint<float, 3>>(false);
	testKdTreeKNearestIndex<TestPoint<double, 3>>(false);
	testKdTreeKNearestIndex<TestPoint<long double, 3>>(false);

    cout << "Test KNearest (from Index) in 3D..." << endl;
	testKdTreeKNearestIndex<TestPoint<float, 4>>(false);
	testKdTreeKNearestIndex<TestPoint<double, 4>>(false);
	testKdTreeKNearestIndex<TestPoint<long double, 4>>(false);
}
