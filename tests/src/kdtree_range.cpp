/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"
#include "../common/kdtree_utils.h"

#include "Ponca/src/SpatialPartitioning/KdTree/kdTreeQuery.h"

using namespace Ponca;

template<typename DataPoint>
void testKdTreeRangeIndex(bool quick = true)
{
	using Scalar = DataPoint::Scalar;
	using VectorContainer = DataPoint::VectorContainer;
	using VectorType = DataPoint::VectorType;

	const int N = quick ? 100 : 100000;
	auto points = VectorContainer(N);
	std::generate(points.begin(), points.end(), []() {return VectorType::Random(); });

	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);

	int seed = 0;

	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		const Scalar r = static_cast<Scalar>((std::rand()) / RAND_MAX * 2.5);
		results.clear();
		KdTreeRangeIndexQuery<DataPoint> rangeIndexQuery = structure.range_neighbors(i, r);
		for (KdTreeRangeIndexIterator<DataPoint> j = rangeIndexQuery.begin(); j != rangeIndexQuery.end(); j++) {
			results.push_back(*j);
		}
		bool res = check_range_neighbors<Scalar, VectorContainer>(points, sampling, i, r, results);
		EXPECT_TRUE(res);
	}

}
template<typename DataPoint>
void testKdTreeRangePoint(bool quick = true)
{
	using Scalar = DataPoint::Scalar;
	using VectorContainer = DataPoint::VectorContainer;
	using VectorType = DataPoint::VectorType;

	const int N = quick ? 100 : 100000;
	auto points = VectorContainer(N);
	std::generate(points.begin(), points.end(), []() {return VectorType::Random(); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

	std::vector<int> results;

	for (int i = 0; i < N; ++i)
	{
		const Scalar r = static_cast<Scalar>((std::rand()) / RAND_MAX * 2.5);
		VectorType point = VectorType::Random();
		results.clear();

		KdTreeRangePointQuery<DataPoint> rangePointQuery = structure.range_neighbors(point, r);
		for (KdTreeRangePointIterator<DataPoint> j = rangePointQuery.begin(); j != rangePointQuery.end(); j++) {
			results.push_back(*j);
		}

		bool res = check_range_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, point, r, results);
		EXPECT_TRUE(res);
	}
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

	testKdTreeRangePoint<TestPoint<float, 3>>(false);
	testKdTreeRangePoint<TestPoint<double, 3>>(false);
	testKdTreeRangePoint<TestPoint<long double, 3>>(false);

	testKdTreeRangePoint<TestPoint<float, 4>>(false);
	testKdTreeRangePoint<TestPoint<double, 4>>(false);
	testKdTreeRangePoint<TestPoint<long double, 4>>(false);

	testKdTreeRangeIndex<TestPoint<float, 3>>(false);
	testKdTreeRangeIndex<TestPoint<double, 3>>(false);
	testKdTreeRangeIndex<TestPoint<long double, 3>>(false);

	testKdTreeRangeIndex<TestPoint<float, 4>>(false);
	testKdTreeRangeIndex<TestPoint<double, 4>>(false);
	testKdTreeRangeIndex<TestPoint<long double, 4>>(false);
}
