/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeNearestPointQuery.h>
#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeKNearestPointQuery.h>
#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeRangePointQuery.h>

#include <Ponca/src/SpatialPartitioning/KdTree/Iterator/KdTreeRangePointIterator.h>
#include <Ponca/src/SpatialPartitioning/KdTree/Iterator/KdTreeKNearestPointIterator.h>

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <chrono>
#include <Eigen/Dense>

//#include <PCP/Math/VectorArray.h>

#include <chrono>
using namespace Eigen;

using namespace std;

using namespace Ponca;

using uint = unsigned int;
using Scalar = float;


//Dans le fichier de test de thibault has_duplicate sert a vérifier les données d'entrée pour empecher les duplicatas
// et retourne un erreur si c'est le cas 

TEST(KdTree_nearest_neighbor_point, sampling)
{
    const int N = quick ? 100 : 10000;
    auto points = std::make_shared<Vector3Array>(N);
    std::generate(points->begin(), points->end(), []() {return Vector3::Random(); });

    std::vector<int> indices(N);
    std::vector<int> sampling(N / 2);
    std::iota(indices.begin(), indices.end(), 0);
    std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(pcptest::seed));

    KdTree structure(points, sampling);

    std::vector<int> results;
    for (int i = 0; i < N; ++i)
    {
        Vector3 point = Vector3::Random();
        results.clear();
        for (int j : structure.nearest_neighbor(point))
        {
            results.push_back(j);
        }
        EXPECT_EQ(int(results.size()), 1);
        EXPECT_TRUE(check_nearest_neighbors(*points, sampling, point, results.front()));
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPositionNormal<Scalar, Dim> Point;

	typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
	typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

	cout << "Testing with parabola..." << endl;
	for (int i = 0; i < g_repeat; ++i)
	{
		CALL_SUBTEST((testFunction<Point, WeightSmoothFunc>()));
		CALL_SUBTEST((testFunction<Point, WeightConstantFunc>()));
	}
	cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

	typedef Eigen::Matrix <Scalar, 3, 1> Vector3;
	const Ponca::KdTree* tree = new Ponca::KdTree();
	//KdTreeKNearestPointQuery<Vector3> q;
	//q.begin();
	KdTreeRangePointIterator<Vector3> t;
}
