/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \file Test general properties of the KdTree
*/

#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../common/has_duplicate.h"
#include "../common/kdtree_utils.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTreeTraits.h>

using namespace Ponca;


template<typename DataPoint>
void testKdtreeWithDuplicate()
{
    using Scalar = typename DataPoint::Scalar;
    using VectorContainer = typename KdTreeSparse<DataPoint>::PointContainer;
    using VectorType = typename DataPoint::VectorType;

    // Number of point samples in each KdTree leaf
#ifdef PONCA_DEBUG
    const int cellSize = 6;
    const int nbCells = 2;
    const int N = nbCells*cellSize;
#else
    const int cellSize = 64;
    const int nbCells = 100;
    const int N = nbCells*cellSize;
#endif

    const Scalar r = 0.001;

    auto test_tree = [r] (const auto& points, const auto&indices, const int cellSize) -> void
    {
        KdTreeSparse<DataPoint> tree(points, indices, cellSize);

#ifndef PONCA_DEBUG
#pragma omp parallel for default(none) shared(tree, points, indices, g_test_stack, r)
#endif
        for (int i = 0; i < points.size(); ++i)
        {
            VectorType point = points[i].pos();//VectorType::Random(); // values between [-1:1]
            std::vector<int> results;

            for (int j : tree.range_neighbors(point, r)) {
                results.push_back(j);
            }

            bool res = check_range_neighbors<Scalar, VectorType, VectorContainer>(points, indices, point, r, results, true);
            VERIFY(res);
        }
    };


    // Generate N random points
    typename KdTreeDense<DataPoint>::IndexContainer ids(N);
    std::iota(ids.begin(), ids.end(), 0);

    auto points = VectorContainer(N);
    std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });

    // Test on 100% random points
    {
        test_tree(points, ids, cellSize);
    }

    // Generate a small part of duplicates by extending the index container
    {
        const int nbDuplicates = N;

        ids.resize(nbDuplicates+N);
        std::generate(ids.begin()+N, ids.end(), [N]() {return Eigen::internal::random<int>(0,N-1); });

        test_tree(points, ids, cellSize);
    }

    // Generate duplicated coordinates samples TODO
    //    {
//        const int nbDuplicates = N/10;
//        const int nbUniques = N;
//
//        auto points = VectorContainer(nbUniques);
//        std::generate(points.begin(), points.end(), []() {return DataPoint(VectorType::Random()); });
//
//        typename KdTreeDense<DataPoint>::IndexContainer ids(nbUniques);
//        std::iota(ids.begin(), ids.end(), 0);
//        ids.resize(nbDuplicates*nbUniques);
//
//        for (int i = 1; i < nbDuplicates; ++i)
//        {
//            std::copy(ids.begin(), ids.begin() + nbUniques, ids.begin()+(nbUniques*i));
//        }
//        test_tree(points, ids, cellSize);
//    }

}

template<typename NodeType>
void testKdTreeNode() {
    std::vector<NodeType> buffer;

    buffer.resize(10);

    // simple predicate that only check if a node is a leaf or not
    auto nodePredicate = [](const NodeType& n1, const NodeType& n2) -> bool {
        return n1.is_leaf() == n2.is_leaf();
    };

    auto checkProperties = [nodePredicate](std::vector<NodeType>& buf, bool targetLeafState) -> void{
        // Check that references works well:
        for (auto& b : buf ) b.set_is_leaf(targetLeafState);
        for (const auto& b : buf ) VERIFY( b.is_leaf() == targetLeafState );

        // Check that copies are working well
        std::vector<NodeType> other;
        other.reserve(buf.size());
        other = buf;
        VERIFY(std::equal(buf.cbegin(), buf.cend(), other.cbegin(), other.cend(), nodePredicate));

        // Check that reallocation works well
        other.resize(buf.size()*2);
        VERIFY(std::equal(buf.cbegin(), buf.cend(), other.cbegin(), other.cbegin()+buf.size(), nodePredicate));
    };

    checkProperties(buffer, true);
    checkProperties(buffer, false);
}

int main(int argc, char** argv)
{
    if (!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    using PointType = TestPoint<float, 3>;
    using KdTreeTraits = KdTreeDefaultTraits<PointType>;

    cout << "Test KdTreeDefaultNode" << endl;
    testKdTreeNode<typename KdTreeTraits::NodeType>();

    cout << "Test KdTreeRange with large number of duplicated points" << endl;
    testKdtreeWithDuplicate<PointType>();

}
