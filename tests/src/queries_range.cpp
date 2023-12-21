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
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>

using namespace Ponca;

template<typename DataPoint, bool SampleKdTree = true>
void testKdTreeRangeIndex(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTree<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 5000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

    /// [KdTree pointer usage]
    // Abstract pointer type that can receive KdTreeSparse or KdTreeDense objects
    KdTree<DataPoint> *kdtree {nullptr};
    /// [KdTree pointer usage]

    std::vector<int> sampling; // we need sampling for GT computation
    if(SampleKdTree){
        std::vector<int> indices(N);
        std::iota(indices.begin(), indices.end(), 0);

        sampling.resize(N / 2);

        int seed = 0;
        std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

        /// [KdTree assign sparse]
        // assign sparse
        kdtree = new KdTreeSparse<DataPoint> (points, sampling);
        /// [KdTree assign sparse]

    } else {
        sampling.resize(N);
        std::iota(sampling.begin(), sampling.end(), 0);
        /// [KdTree assign dense]
        // assign dense
        kdtree = new KdTreeDense<DataPoint> (points);
        /// [KdTree assign dense]
    }

#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        Scalar r = Eigen::internal::random<Scalar>(0., 0.5);
        std::vector<int> resultsTree;

        for (int j : kdtree->range_neighbors(i, r)) {
            resultsTree.push_back(j);
        }
        if( SampleKdTree ) {
            bool resTree = check_range_neighbors<Scalar, VectorContainer>(points, sampling, i, r, resultsTree);
            VERIFY(resTree);
        }
        else {
            bool resTree = check_range_neighbors<Scalar, VectorContainer>(points, sampling, i, r, resultsTree);
            VERIFY(resTree);
        }
    }

    if( ! SampleKdTree ){
        Ponca::KnnGraph<DataPoint> knnGraph(*kdtree, N/4); // we need a large graph, otherwise we might miss some points
                                                           // (which is the goal of the graph: to replace full euclidean
                                                           // collection by geodesic-like region growing bounded by
                                                           // the euclidean ball).
#pragma omp parallel for
        for (int i = 0; i < N; ++i)
        {
            Scalar r = Eigen::internal::random<Scalar>(0., 0.5);
            std::vector<int> resultsGraph;

            for (int j : knnGraph.range_neighbors(i, r)) {
                resultsGraph.push_back(j);
            }
            bool resGraph = check_range_neighbors<Scalar, VectorContainer>(points, sampling, i, r, resultsGraph);
            VERIFY(resGraph);
        }
    }

    delete kdtree;
}
template<typename DataPoint>
void testKdTreeRangePoint(bool quick = true)
{
	using Scalar = typename DataPoint::Scalar;
	using VectorContainer = typename KdTreeSparse<DataPoint>::PointContainer;
	using VectorType = typename DataPoint::VectorType;

	const int N = quick ? 100 : 10000;
	auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

	int seed = 0;

    /// [Kdtree sampling construction]
	std::vector<int> indices(N);
	std::vector<int> sampling(N / 2);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));

	KdTreeSparse<DataPoint> structure(points, sampling);
    /// [Kdtree sampling construction]

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

#ifndef NDEBUG
    bool quick = true;
#else
    bool quick = false;
#endif

    cout << "Test KdTreeRange (from Point) in 3D..." << endl;
	testKdTreeRangePoint<TestPoint<float, 3>>(quick);
	testKdTreeRangePoint<TestPoint<double, 3>>(quick);
	testKdTreeRangePoint<TestPoint<long double, 3>>(quick);

    cout << "Test KdTreeRange (from Point) in 4D..." << endl;
	testKdTreeRangePoint<TestPoint<float, 4>>(quick);
	testKdTreeRangePoint<TestPoint<double, 4>>(quick);
	testKdTreeRangePoint<TestPoint<long double, 4>>(quick);

    cout << "Test Range Queries (from Index) using KnnGraph and Kdtree in 3D... (without subsampling)" << endl;
    testKdTreeRangeIndex<TestPoint<float, 3>, false>(quick);
    testKdTreeRangeIndex<TestPoint<double, 3>, false>(quick);
    testKdTreeRangeIndex<TestPoint<long double, 3>, false>(quick);

    cout << "Test KdTreeRange (from Index) in 3D... (with subsampling)" << endl;
    testKdTreeRangeIndex<TestPoint<float, 3>>(quick);
    testKdTreeRangeIndex<TestPoint<double, 3>>(quick);
    testKdTreeRangeIndex<TestPoint<long double, 3>>(quick);

    cout << "Test KdTreeRange (from Index) in 4D..." << endl;
    testKdTreeRangeIndex<TestPoint<float, 4>>(quick);
    testKdTreeRangeIndex<TestPoint<double, 4>>(quick);
    testKdTreeRangeIndex<TestPoint<long double, 4>>(quick);

    cout << "Test Range Queries (from Index) using KnnGraph and Kdtree in 4D... (without subsampling)" << endl;
    testKdTreeRangeIndex<TestPoint<float, 4>, false>(quick);
    testKdTreeRangeIndex<TestPoint<double, 4>, false>(quick);
    testKdTreeRangeIndex<TestPoint<long double, 4>, false>(quick);
}
