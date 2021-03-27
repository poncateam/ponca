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
void testKdTreeNearestIndex(bool quick = true)
{
	using Scalar = DataPoint::Scalar;
	using VectorContainer = DataPoint::VectorContainer;
	using VectorType = DataPoint::VectorType;

	const int N = quick ? 100 : 100000;
	auto points = VectorContainer(N);
	std::generate(points.begin(), points.end(), []() {return VectorType::Random(); });

	KdTree<DataPoint> structure(points);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		results.clear();
		for (int j : structure.nearest_neighbor(i))
		{
			results.push_back(j);
		}
		EXPECT_EQ(int(results.size()), 1);
		bool res = check_nearest_neighbor<Scalar, VectorContainer>(points, i, results.front());
		EXPECT_TRUE(res);
	}
}

template<typename DataPoint>
void testKdTreeNearestPoint(bool quick = true)
{
	using Scalar = DataPoint::Scalar;
	using VectorContainer = DataPoint::VectorContainer;
	using VectorType = DataPoint::VectorType;

	const int N = quick ? 100 : 100000;
	auto points = VectorContainer(N);
	std::generate(points.begin(), points.end(), []() {return VectorType::Random(); });

	KdTree<DataPoint> structure(points);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		VectorType point = VectorType::Random();
		results.clear();
		for (int j : structure.nearest_neighbor(point))
		{
			results.push_back(j);
		}
		EXPECT_EQ(int(results.size()), 1);
		bool res = check_nearest_neighbor<Scalar, VectorType, VectorContainer>(points, point, results.front());
		EXPECT_TRUE(res);
	}
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

	testKdTreeNearestPoint<TestPoint<float, 3>>(false);
	testKdTreeNearestPoint<TestPoint<double, 3>>(false);
	testKdTreeNearestPoint<TestPoint<long double, 3>>(false);

	testKdTreeNearestPoint<TestPoint<float, 4>>(false);
	testKdTreeNearestPoint<TestPoint<double, 4>>(false);
	testKdTreeNearestPoint<TestPoint<long double, 4>>(false);

	testKdTreeNearestIndex<TestPoint<float, 3>>(false);
	testKdTreeNearestIndex<TestPoint<double, 3>>(false);
	testKdTreeNearestIndex<TestPoint<long double, 3>>(false);

	testKdTreeNearestIndex<TestPoint<float, 4>>(false);
	testKdTreeNearestIndex<TestPoint<double, 4>>(false);
	testKdTreeNearestIndex<TestPoint<long double, 4>>(false);
}
