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


template<typename DataPoint,
		typename VectorContainer,
		typename Scalar = typename DataPoint::Scalar,
		typename VectorType = typename DataPoint::VectorType>
void check_support(
	KdTreeRangePointQuery<DataPoint>& supportQuery, const VectorContainer& points, 
	const vector<VectorType>& queryPoints,
	const std::vector<int>& sampling, Scalar r)
{
	std::vector<int> results;
	for(int q : supportQuery)
	{
		results.push_back(q);
	}
	VectorType center = VectorType::Zero();
	Scalar support = 0;

	for (int i = 0; i < queryPoints.size(); i++)
	{
		center =center + queryPoints[i];
	}

	center /= queryPoints.size();


	for (int i = 0; i < queryPoints.size(); i++)
	{
		Scalar temp = std::sqrt((center - queryPoints[i]).squaredNorm());
		if (support < temp)
		{
			support = temp;
		}
	}

	bool res = check_range_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, center, support + r, results);
	VERIFY(res);
}

template<typename DataPoint>
void testKdTreeSupportRangePoint(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ?  100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);


    //Initialize Queries
    const int queryCount = quick ? 100 : 1000;
    Scalar r = Eigen::internal::random<Scalar>(0.0, 0.01);
    vector<VectorType> queryPoints;

    for(int i = 0;i < queryCount;i++)
    {
		queryPoints.push_back(VectorType::Random()); // values between [-1:1]
    }

    auto supportQuery = structure.support_range_neighbors(queryPoints,r);

	
	check_support<DataPoint,VectorContainer>(supportQuery,points,queryPoints, sampling, r);
}


template<typename DataPoint>
void testKdTreeSupportRangeIndex(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ?  100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	int seed = 0;
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTree<DataPoint> structure(points, sampling);


    //Initialize Queries
    const int queryCount = quick ? 100 : 1000;
    Scalar r = Eigen::internal::random<Scalar>(0.0, 0.01);
    vector<int> queryPoints;

    for(int i = 0;i < queryCount;i++)
    {
		queryPoints.push_back(rand()%N); // values between [-1:1]
    }

    auto supportQuery = structure.support_range_neighbors(queryPoints,r);

	
	check_support<DataPoint,VectorContainer>(supportQuery,points,queryPoints, sampling, r);
}

template<typename DataPoint,
		typename VectorContainer,
		typename Scalar = typename DataPoint::Scalar,
		typename VectorType = typename DataPoint::VectorType>
void check_support(
	KdTreeRangePointQuery<DataPoint>& supportQuery, const VectorContainer& points, 
	const vector<int>& queryPoints,
	const std::vector<int>& sampling, Scalar r)
{
	std::vector<int> results;
	for(int q : supportQuery)
	{
		results.push_back(q);
	}
	VectorType center = VectorType::Zero();
	Scalar support = 0;

	for (int i = 0; i < queryPoints.size(); i++)
	{
		center =center + points[queryPoints[i]].pos();
	}

	center /= queryPoints.size();


	for (int i = 0; i < queryPoints.size(); i++)
	{
		Scalar temp = std::sqrt((center - points[queryPoints[i]].pos()).squaredNorm());
		if (support < temp)
		{
			support = temp;
		}
	}

	bool res = check_range_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, center, support + r, results);
	VERIFY(res);
}

int main(int argc, char** argv)
{
	if (!init_testing(argc, argv))
	{
		return EXIT_FAILURE;
	}

    cout << "Test Stacked Range (from Point) in 3D..." << endl;
	testKdTreeSupportRangePoint<TestPoint<float, 3>>(true);
	testKdTreeSupportRangePoint<TestPoint<double, 3>>(true);
	testKdTreeSupportRangePoint<TestPoint<long double, 3>>(true);

    cout << "Test Stacked Range (from Point) in 4D..." << endl;
	testKdTreeSupportRangePoint<TestPoint<float, 4>>(true);
	testKdTreeSupportRangePoint<TestPoint<double, 4>>(true);
	testKdTreeSupportRangePoint<TestPoint<long double, 4>>(true);

    cout << "Test Stacked Range (from Index) in 3D..." << endl;
	testKdTreeSupportRangeIndex<TestPoint<float, 3>>(true);
	testKdTreeSupportRangeIndex<TestPoint<double, 3>>(true);
	testKdTreeSupportRangeIndex<TestPoint<long double, 3>>(true);

    cout << "Test Stacked Range (from Index) in 4D..." << endl;
	testKdTreeSupportRangeIndex<TestPoint<float, 4>>(true);
	testKdTreeSupportRangeIndex<TestPoint<double, 4>>(true);
	testKdTreeSupportRangeIndex<TestPoint<long double, 4>>(true);
}
