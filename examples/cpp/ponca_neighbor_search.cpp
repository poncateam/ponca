/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <random>
#include <Ponca/SpatialPartitioning>
#include <Eigen/Core>

struct DataPoint
{
    enum {Dim = 3};
    using Scalar = float;
    using VectorType = Eigen::Vector<Scalar,Dim>;
    [[nodiscard]] inline const auto& pos() const {return m_pos;}
    VectorType m_pos;
};

int main()
{
    // generate N random points
    constexpr int N {100000};


    //////////////////////////////////////////////////////////////
    ////////////////////// Tree construction /////////////////////
    //////////////////////////////////////////////////////////////
    /// [Kdtree Dense construction]
    std::vector<DataPoint> points(N);
    std::generate(points.begin(), points.end(), [](){
        return DataPoint{100 * DataPoint::VectorType::Random()};
    });
    // build the k-d tree
    Ponca::KdTreeDense<DataPoint> kdtreeDense(points);
    /// [Kdtree Dense construction]

    int seed = 0;
    /// [Kdtree Sampling construction]
    std::vector<int> indices(N);
    std::vector<int> sampling(N / 2);
    std::iota(indices.begin(), indices.end(), 0);
    std::sample(indices.begin(), indices.end(), sampling.begin(), N / 2, std::mt19937(seed));
    Ponca::KdTreeSparse<DataPoint> kdtreeSparse(points, sampling);
    /// [Kdtree Sampling construction]

    /// [KdTree pointer usage]
    // Abstract pointer type that can receive KdTreeSparse or KdTreeDense objects
    Ponca::KdTree<DataPoint> *kdtree {nullptr};
    /// [KdTree pointer usage]

    /// [KdTree assign sparse]
    // Assign sparse
    kdtree = &kdtreeSparse; // Or kdtree = new Ponca::KdTreeSparse<DataPoint> (points, sampling);
    /// [KdTree assign sparse]

    /// [KdTree assign dense]
    // Assign dense
    kdtree = &kdtreeDense; // Or kdtree = new Ponca::KdTreeDense<DataPoint> (points, sampling);
    /// [KdTree assign dense]

    // neighbor searches are done below from these arbitrary index and point
    const int query_idx = 10;
    const DataPoint::VectorType query_pt{-10.0, 0.5, 75.0};

    const int k = 10;
    /// [KnnGraph construction]
    Ponca::KnnGraph<DataPoint> knnGraph(kdtreeDense, k);
    /// [KnnGraph construction]

    //////////////////////////////////////////////////////////////
    ////////////////// Nearest neighbors search //////////////////
    //////////////////////////////////////////////////////////////
    std::cout << "The nearest neighbor of the point at index " << query_idx << " is at index "
              << *kdtree->nearest_neighbor(query_idx).begin() << std::endl;
    std::cout << "The nearest neighbor of the point (" << query_pt.transpose() << ") is at index "
              << *kdtree->nearest_neighbor(query_pt).begin() << std::endl;

    /// [Kdtree k-nearest neighbor index search]
    std::cout << "The " << k << "-nearest neighbors of the point at index " << query_idx << " are at indices: ";
    for(int neighbor_idx : kdtreeDense.k_nearest_neighbors(query_idx, k)) { // Iterates over the neighbors of query_idx
        std::cout << neighbor_idx << ", ";
    }
    /// [Kdtree k-nearest neighbor index search]
    std::cout << std::endl;

    /// [Kdtree k-nearest neighbor position search]
    std::cout << "The " << k << "-nearest neighbors of the point (" << query_pt.transpose() << ") are at indices: ";
    for(auto neighbor_idx : kdtreeDense.k_nearest_neighbors(query_pt, k)) {
        std::cout << neighbor_idx << ", ";
    }
    /// [Kdtree k-nearest neighbor position search]
    std::cout << std::endl;

    /// [KnnGraph k-nearest neighbor index search]
    std::cout << "The nearest neighbors of the point at index " << query_idx << " are at indices: ";
    for(int neighbor_idx : knnGraph.k_nearest_neighbors(query_idx)) {
        std::cout << neighbor_idx << ", ";
    }
    /// [KnnGraph k-nearest neighbor index search]
    std::cout << std::endl;

    //////////////////////////////////////////////////////////////
    /////////////////// Range neighbors search ///////////////////
    //////////////////////////////////////////////////////////////
    const int second_query_idx = 5;
    /// [Kdtree range neighbors index mutable search]
    constexpr DataPoint::Scalar radius = 5.25;
    auto rangeNeighbors = kdtreeDense.range_neighbors(query_idx, radius); // First query

    std::cout << "The neighbors of the point at index " << second_query_idx << " at a distance " << radius << " are at indices: ";
    for(int neighbor_idx : rangeNeighbors(second_query_idx, radius)) { // Iterate over a new range_neighbors query
        std::cout << neighbor_idx << ", ";
    }
    /// [Kdtree range neighbors index mutable search]
    std::cout << std::endl;
    /// [Kdtree range neighbors position mutable search]
    auto posRangeNeighbors = kdtreeDense.range_neighbors_empty_position(); // Empty query

    std::cout << "The neighbors of the point (" << query_pt.transpose() << ") at a distance " << radius << " are at indices: ";
    for(auto neighbor_idx : posRangeNeighbors(query_pt, radius)) {
        std::cout << neighbor_idx << ", ";
    }
    /// [Kdtree range neighbors position mutable search]
    std::cout << std::endl;

    /// [KnnGraph range neighbors index mutable search]
    auto knnRangeNeighbors = knnGraph.range_neighbors_empty_index(); // Empty query

    std::cout << "The neighbors of the point (" << query_pt.transpose() << ") at a distance " << radius << " are at indices: ";
    for(auto neighbor_idx : knnRangeNeighbors(query_idx, radius)) { // Iterates over the neighbors of query_idx
        std::cout << neighbor_idx << ", ";
    }
    /// [KnnGraph range neighbors index mutable search]
    std::cout << std::endl;

    return 1;
}
