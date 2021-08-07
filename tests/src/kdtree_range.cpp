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
void testKdTreeRangeIndex(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);

	int seed = 0;

	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
        Scalar r = Eigen::internal::random<Scalar>(0., 0.5);
        std::vector<int> results;

		for (int j : structure.range_neighbors(i, r)) {
			results.push_back(j);
		}
		bool res = check_range_neighbors<Scalar, VectorContainer>(points, sampling, i, r, results);
        VERIFY(res);
	}

}
template<typename DataPoint>
void testKdTreeRangePoint(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
        Scalar r = Eigen::internal::random<Scalar>(0., 0.5);
		VectorType point = VectorType::Random(); // values between [-1:1]
        std::vector<int> results;

		for (int j : structure.range_neighbors(point, r)) {
			results.push_back(j);
		}

		bool res = check_range_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, point, r, results);
		VERIFY(res);
	}
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

    cout << "Test KdTreeRange (from Point) in 3D..." << endl;
	testKdTreeRangePoint<TestPoint<float, 3>>(false);
	testKdTreeRangePoint<TestPoint<double, 3>>(false);
	testKdTreeRangePoint<TestPoint<long double, 3>>(false);

    cout << "Test KdTreeRange (from Point) in 4D..." << endl;
	testKdTreeRangePoint<TestPoint<float, 4>>(false);
	testKdTreeRangePoint<TestPoint<double, 4>>(false);
	testKdTreeRangePoint<TestPoint<long double, 4>>(false);

    cout << "Test KdTreeRange (from Index) in 3D..." << endl;
	testKdTreeRangeIndex<TestPoint<float, 3>>(false);
	testKdTreeRangeIndex<TestPoint<double, 3>>(false);
	testKdTreeRangeIndex<TestPoint<long double, 3>>(false);

    cout << "Test KdTreeRange (from Index) in 4D..." << endl;
	testKdTreeRangeIndex<TestPoint<float, 4>>(false);
	testKdTreeRangeIndex<TestPoint<double, 4>>(false);
	testKdTreeRangeIndex<TestPoint<long double, 4>>(false);
}
