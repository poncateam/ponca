/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"

#include "Ponca/src/SpatialPartitioning/KdTree/kdTreeQuery.h"

#include <chrono>
#include <Eigen/Dense>
#include <cassert>

#include <chrono>
#include <random>
#include <iterator>
#include <algorithm>


#define EXPECT_EQ(a,b) assert(a==b)
#define EXPECT_TRUE(a) assert(a==True)

using namespace Eigen;
//using namespace Ponca;
using namespace std;



using uint = unsigned int;

using Vec3 = Eigen::Matrix<float, 3, 1>;
using Vector3Array = std::vector<Vec3>;




//Dans le fichier de test de thibault has_duplicate sert a vérifier les données d'entrée pour empecher les duplicatas
// et retourne un erreur si c'est le cas 


template <typename Scalar, typename Vector3>
class _Point {
public:
	typedef Scalar Scalar;
	typedef Vector3 VectorType;
};




bool check_k_nearest_neighbors(const Vector3Array& points, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	float max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if (idx == index) continue;

		float dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}

bool check_nearest_neighbors(const Vector3Array& points, int index, int nearest)
{
	return check_k_nearest_neighbors(points, index, 1, { nearest });
}

bool check_k_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	float max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		if (idx == index) continue;

		float dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}

bool check_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, int index, int nearest)
{
	return check_k_nearest_neighbors(points, sampling, index, 1, { nearest });
}


template<typename Vector3>
void callSubTests()
{
	const int N = 10000;// quick ? 100 : 10000
	auto points = std::make_shared<Vector3Array>(N);
	//std::generate(points->begin(), points->end(), []() {return Vector3::Random(); });


	//std::vector<int> indices(N);
	//std::vector<int> sampling(N / 2);
	//std::iota(indices.begin(), indices.end(), 0);
	
	//NEED C++17
	//std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(pcptest::seed));

	using Point = _Point<float, Vector3Array>;


	Ponca::KdTree<Point>* structure = new Ponca::KdTree<Point>(points);

	//std::vector<int> results;
	//for (int i = 0; i < N; ++i)
	//{
	//	const float r = float(std::rand()) / RAND_MAX * 2.5;
	//	results.clear();
	//	for (int j : structure.range_neighbors(i, r))
	//	{
	//		results.push_back(j);
	//	}
	//	//EXPECT_EQ(int(results.size()), 1);
	//	//EXPECT_TRUE(check_nearest_neighbors(*points, i, results.front()));
	//}
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