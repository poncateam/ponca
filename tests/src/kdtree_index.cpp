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
	const int N = 10;// quick ? 100 : 10000
	auto points = std::make_shared<Vector3Array>(N);
	std::generate(points->begin(), points->end(), []() {return Vector3::Random(); });


	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	
	int seed = 0;

	//NEED C++17
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));
	
	using Point = _Point<float, Vector3Array>;

	KdTree<Point> structure(points, sampling);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		const float r = float(std::rand()) / RAND_MAX * 2.5;
		results.clear();
		for (int j : structure.range_neighbors(i, r))
		{
			results.push_back(j);
		}
		//EXPECT_EQ(int(results.size()), 1);
		//EXPECT_TRUE();
		EXPECT_TRUE(check_range_neighbors(*points,sampling, i, r,results));
	}
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

	//typedef Eigen::Matrix <Scalar, 3, 1> Vector3;
	//const Ponca::KdTree* tree = new Ponca::KdTree();
	//KdTreeNearestIndexQuery q;
	//q.begin();

	callSubTests<Vec3>();

}