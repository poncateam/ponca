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
void testKdTreeRangeIndex(bool quick)
{
	using S = DataPoint::Scalar;
	using V = DataPoint::Vector;
	using VT = DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = std::make_shared<V>(N);
	std::generate(points->begin(), points->end(), []() {return VT::Random(); });

	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);

	int seed = 0;

	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		const S r = static_cast<S>((std::rand()) / RAND_MAX * 2.5);
		results.clear();
		KdTreeRangeIndexQuery<DataPoint> rangeIndexQuery = structure.range_neighbors(i, r);
		for (KdTreeRangeIndexIterator<DataPoint> j = rangeIndexQuery.begin(); j != rangeIndexQuery.end(); j++) {
			results.push_back(*j);
		}
		bool res = check_range_neighbors<S, V>(*points, sampling, i, r, results);
		EXPECT_TRUE(res);
	}

}
template<typename DataPoint>
void testKdTreeRangePoint(bool quick)
{
	using S = DataPoint::Scalar;
	using V = DataPoint::Vector;
	using VT = DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = std::make_shared<V>(N);
	std::generate(points->begin(), points->end(), []() {return VT::Random(); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);

	std::vector<int> results;

	for (int i = 0; i < N; ++i)
	{
		const S r = static_cast<S>((std::rand()) / RAND_MAX * 2.5);
		VT point = VT::Random();
		results.clear();

		KdTreeRangePointQuery<DataPoint> rangePointQuery = structure.range_neighbors(point, r);
		for (KdTreeRangePointIterator<DataPoint> j = rangePointQuery.begin(); j != rangePointQuery.end(); j++) {
			results.push_back(*j);
		}

		bool res = check_range_neighbors<S, VT, V>(*points, sampling, point, r, results);
		EXPECT_TRUE(res);
	}
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

	testKdTreeRangePoint<_Point<float, 3>>(false);
	testKdTreeRangePoint<_Point<double, 3>>(false);
	testKdTreeRangePoint<_Point<long double, 3>>(false);

	testKdTreeRangePoint<_Point<float, 4>>(false);
	testKdTreeRangePoint<_Point<double, 4>>(false);
	testKdTreeRangePoint<_Point<long double, 4>>(false);

	testKdTreeRangeIndex<_Point<float, 3>>(false);
	testKdTreeRangeIndex<_Point<double, 3>>(false);
	testKdTreeRangeIndex<_Point<long double, 3>>(false);

	testKdTreeRangeIndex<_Point<float, 4>>(false);
	testKdTreeRangeIndex<_Point<double, 4>>(false);
	testKdTreeRangeIndex<_Point<long double, 4>>(false);
}
