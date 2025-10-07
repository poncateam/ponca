/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <Ponca/SpatialPartitioning>
#include <Eigen/Core>

struct DataPoint
{
    enum {Dim = 3};
    using Scalar = float;
    using VectorType = Eigen::Vector<Scalar,Dim>;
    inline const auto& pos() const {return m_pos;}
    VectorType m_pos;
};

int main()
{
    // generate N random points
    constexpr int N {100000};
    std::vector<DataPoint> points(N);
    std::generate(points.begin(), points.end(), [](){
        return DataPoint{100 * DataPoint::VectorType::Random()};});

    // build the k-d tree
    const Ponca::KdTreeDense<DataPoint> kdtree(points);

    // neighbor searches are done below from these arbitrary index and point
    const int query_idx = 10;
    const DataPoint::VectorType query_pt{-10.0, 0.5, 75.0};

    //
    // nearest neighbor search
    //
    std::cout << "the nearest neighbor of the point at index " << query_idx << " is at index "
              << *kdtree.nearest_neighbor(query_idx).begin() << std::endl;
    std::cout << "the nearest neighbor of the point (" << query_pt.transpose() << ") is at index "
              << *kdtree.nearest_neighbor(query_pt).begin() << std::endl;

    //
    // k-nearest neighbor search
    //
    constexpr int k = 10;
    std::cout << "the " << k << "-nearest neighbors of the point at index " << query_idx << " are at indices: ";
    for(int neighbor_idx : kdtree.k_nearest_neighbors(query_idx, k)) {
        std::cout << neighbor_idx << ", ";
    }
    std::cout << std::endl;
    std::cout << "the " << k << "-nearest neighbors of the point (" << query_pt.transpose() << ") are at indices: ";
    for(auto neighbor_idx : kdtree.k_nearest_neighbors(query_pt, k)) {
        std::cout << neighbor_idx << ", ";
    }
    std::cout << std::endl;

    //
    // range neighbor search
    //
    constexpr DataPoint::Scalar radius = 5.25;
    std::cout << "the neighbors of the point at index " << query_idx << " at a distance " << radius << " are at indices: ";
    for(int neighbor_idx : kdtree.range_neighbors(query_idx, radius)) {
        std::cout << neighbor_idx << ", ";
    }
    std::cout << std::endl;
    std::cout << "the neighbors of the point (" << query_pt.transpose() << ") at a distance " << radius << " are at indices: ";
    for(auto neighbor_idx : kdtree.range_neighbors(query_pt, radius)) {
        std::cout << neighbor_idx << ", ";
    }
    std::cout << std::endl;

    return 1;
}
