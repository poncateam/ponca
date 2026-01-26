/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <optional>
#include <Ponca/SpatialPartitioning>
#include <Ponca/src/Common/pointTypes.h>

#include <Eigen/Core>

using DataPoint = Ponca::PointPosition<float, 3>;

//! [CustomInnerNodeDefinition]
template <typename NodeIndex, typename Scalar, int DIM, typename _AabbType = Eigen::AlignedBox<Scalar, DIM>>
struct MyKdTreeInnerNode : public Ponca::KdTreeDefaultInnerNode<NodeIndex, Scalar, DIM> {
    using AabbType = _AabbType;
    AabbType m_aabb{};
};
//! [CustomInnerNodeDefinition]

//! [CustomNodeDefinition]
template <typename Index, typename NodeIndex, typename DataPoint, typename LeafSize = Index>
struct MyKdTreeNode : Ponca::KdTreeCustomizableNode<Index, NodeIndex, DataPoint, LeafSize,
        MyKdTreeInnerNode<NodeIndex, typename DataPoint::Scalar, DataPoint::Dim>> {

    using Base = Ponca::KdTreeCustomizableNode<Index, NodeIndex, DataPoint, LeafSize,
            MyKdTreeInnerNode<NodeIndex, typename DataPoint::Scalar, DataPoint::Dim>>;
    using AabbType  = typename Base::AabbType;

    void configure_range(Index start, Index size, const AabbType &aabb)
    {
        Base::configure_range(start, size, aabb);
        if (! Base::is_leaf() )
        {
            Base::getAsInner().m_aabb = aabb;
        }
    }
    [[nodiscard]] inline std::optional<AabbType> getAabb() const {
        if (! Base::is_leaf())
            return Base::getAsInner().m_aabb;
        else
            return std::optional<AabbType>();
    }
};
//! [CustomNodeDefinition]

int main()
{
    // generate N random points
    constexpr int N {100000};
    std::vector<DataPoint> points(N);
    std::generate(points.begin(), points.end(), [](){
        return DataPoint{100 * DataPoint::VectorType::Random()};});

//! [KdTreeTypeWithCustomNode]
    using CustomKdTree = Ponca::KdTreeDenseBase<Ponca::KdTreeDefaultTraits<DataPoint,MyKdTreeNode>>;
//! [KdTreeTypeWithCustomNode]

    // build the k-d tree
    const CustomKdTree kdtree(points);

    // neighbor searches are done below from these arbitrary index and point
    const int query_idx = 10;
    const DataPoint::VectorType query_pt{-10.0, 0.5, 75.0};

    //
    // nearest neighbor search
    //
    std::cout << "the nearest neighbor of the point at index " << query_idx << " is at index "
              << *kdtree.nearestNeighbor(query_idx).begin() << std::endl;
    std::cout << "the nearest neighbor of the point (" << query_pt.transpose() << ") is at index "
              << *kdtree.nearestNeighbor(query_pt).begin() << std::endl;

    //! [ReadCustomProperties]
    auto bbox = kdtree.nodes()[0].getAabb();
    if (bbox) {
        std::cout << "Root bounding box is as follows: \n"
                  << "  Center:   " << bbox->center()
                  << "  Diagonal: " << bbox->diagonal()
                  << std::endl;
    }
    //! [ReadCustomProperties]

    return 1;
}
