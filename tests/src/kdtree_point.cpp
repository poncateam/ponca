/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"
#include "../common/kdtree_utils.h"

#include "Ponca/src/SpatialPartitioning/KdTree/kdTreeQuery.h"

using namespace Ponca;

template<typename Vector3>
void callSubTests()
{
	const int N = 10;//quick ? 100 : 10000;
	auto points = std::make_shared<Vector3Array>(N);
	std::generate(points->begin(), points->end(), []() {return Vector3::Random(); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	using Point = _Point<float, Vector3Array>;
	KdTree<Point> structure(points, sampling);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		const float r = float(std::rand()) / RAND_MAX * 2.5;
		Vector3 point = Vector3::Random();
		results.clear();


		KdTreeRangePointQuery<Point> rangePointQuery = structure.range_neighbors(point, r);
		for (KdTreeRangePointIterator<Point> j = rangePointQuery.begin(); j != rangePointQuery.end(); j++) {
			results.push_back(*j);
		}
/*
		for (int j : structure.range_neighbors(point, r))
		{
			results.push_back(j);
		}*/
		EXPECT_TRUE(check_range_neighbors(*points, sampling, point, r, results));
	}
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }
	callSubTests<Vec3>();
}