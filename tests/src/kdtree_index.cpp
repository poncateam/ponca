/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"

#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeNearestIndexQuery.h>

#include <chrono>
#include <Eigen/Dense>

//#include <PCP/Math/VectorArray.h>

#include <chrono>
using namespace Eigen;

using namespace std;

using namespace Ponca;

using uint = unsigned int;
using Scalar = float;

using Vec3 = Eigen::Matrix<SPScalar, 3, 1>;
using Vector3Array = std::vector<Vec3>;

//Dans le fichier de test de thibault has_duplicate sert a vérifier les données d'entrée pour empecher les duplicatas
// et retourne un erreur si c'est le cas 


bool check_k_nearest_neighbors(const Vector3Array& points, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		//PCP_DEBUG_ERROR;;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//PCP_DEBUG_ERROR;;
		return false;
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		//PCP_DEBUG_ERROR;;
		return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if (idx == index) continue;

		Scalar dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//PCP_DEBUG_ERROR;;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//PCP_DEBUG_ERROR;;
			return false;
		}
	}
	return true;
}

template<typename Vector3Array>
bool check_k_nearest_neighbors(const Vector3Array& points, int index, int k, const std::vector<int>& neighbors)
{
	return check_k_nearest_neighbors(points, index, 1, { nearest });


}

template<typename Vector3, int N>
void callSubTests()
{
	const int N = quick ? 100 : 10000;
	auto points = std::make_shared<Vector3Array>(N);
	std::generate(points->begin(), points->end(), []() {return Vector3::Random(); });

	KdTree structure(points);

	std::vector<int> results;
	for (int i = 0; i < N; ++i)
	{
		results.clear();
		for (int j : structure.nearest_neighbor(i))
		{
			results.push_back(j);
		}
		EXPECT_EQ(int(results.size()), 1);
		EXPECT_TRUE(check_nearest_neighbors(*points, i, results.front()));
	}
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

	/*typedef Eigen::Matrix <Scalar, 3, 1> Vector3;
	const Ponca::KdTree* tree = new Ponca::KdTree();
	KdTreeNearestIndexQuery q;
	q.begin();*/

	callSubTests<Vec3,20>();

}
