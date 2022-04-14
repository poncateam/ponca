/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./kdTreeNode.h"

#include <Eigen/Eigen>
#include <Eigen/Geometry> // aabb

#include <memory>
#include <vector>
#include <numeric>

#include "../../Common/Assert.h"

#include "Query/kdTreeNearestPointQuery.h"
#include "Query/kdTreeNearestIndexQuery.h"
#include "Query/kdTreeKNearestPointQuery.h"
#include "Query/kdTreeKNearestIndexQuery.h"
#include "Query/kdTreeRangeIndexQuery.h"
#include "Query/kdTreeRangePointQuery.h"

#define PCA_KDTREE_MAX_DEPTH 32

namespace Ponca {

/// \ingroup spatialpartitioning
///
/// \tparam DataPoint
///
/// \todo Better handle sampling: do not store non-selected points (requires to store original indices
template<class _DataPoint>
class KdTree
{
public:
    typedef          _DataPoint            DataPoint;
	typedef typename DataPoint::Scalar     Scalar;  // Scalar given by user
	typedef typename DataPoint::VectorType VectorType;  // VectorType given by user

	typedef typename Eigen::AlignedBox<Scalar, DataPoint::Dim> Aabb; // Intersections

    typedef typename std::vector<DataPoint> PointContainer; // Container for VectorType used inside the KdTree
    typedef typename std::vector<int> IndexContainer; // Container for indices used inside the KdTree
    typedef typename std::vector<KdTreeNode<Scalar>> NodeContainer;  // Container for nodes used inside the KdTree

    inline KdTree():
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
    };

    template<typename PointUserContainer>
    inline KdTree(const PointUserContainer& points): // PointUserContainer => Given by user, transformed to PointContainer
        m_points(PointContainer()),
        m_nodes(NodeContainer()),
        m_indices(IndexContainer()),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        this->build(points);
    };

    template<typename PointUserContainer, typename IndexUserContainer>
    inline KdTree(const PointUserContainer& points, const IndexUserContainer& sampling): // PointUserContainer => Given by user, transformed to PointContainer
                                                                                         // IndexUserContainer => Given by user, transformed to IndexContainer
        m_points(),
        m_nodes(),
        m_indices(),
        m_min_cell_size(64),
        m_leaf_count(0)
    {
        buildWithSampling(points, sampling);
    };

    inline void clear();

    struct DefaultConverter{
        template <typename Input>
        inline void operator()( const Input& i, PointContainer & o ) {
            if constexpr ( std::is_same<Input, PointContainer>::value )
                o = i;
            else
                std::transform(i.cbegin(), i.cend(), std::back_inserter(o),
                               [](const typename Input::value_type &p) -> DataPoint { return DataPoint(p); });
        }
    };

    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \param points Input points
    template<typename PointUserContainer>
    inline void build(const PointUserContainer& points) { build(points, DefaultConverter()); }
    ///
    /// \tparam PointUserContainer Input point container, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling container, transformed to IndexContainer
    /// \param points Input points
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename Converter>
    inline void build(const PointUserContainer& points, Converter c);

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    template<typename PointUserContainer, typename IndexUserContainer>
    inline void buildWithSampling(const PointUserContainer& points,
                                  const IndexUserContainer& sampling)
                                  { buildWithSampling(points, sampling, DefaultConverter());}

    /// \tparam PointUserContainer Input point, transformed to PointContainer
    /// \tparam IndexUserContainer Input sampling, transformed to IndexContainer
    /// \tparam Converter
    /// \param points Input points
    /// \param sampling Indices of points used in the tree
    /// \param c Cast/Convert input point type to DataType
    template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
    inline void buildWithSampling(const PointUserContainer& points,
                                  const IndexUserContainer& sampling,
                                  Converter c);


    template<typename IndexUserContainer>
    inline void rebuild(const IndexUserContainer& sampling); // IndexUserContainer => Given by user, transformed to IndexContainer


    inline bool valid() const;
    inline std::string to_string() const;

    // Accessors ---------------------------------------------------------------
public:
    inline int node_count() const;
    inline int index_count() const;
    inline int point_count() const;
    inline int leaf_count() const;

    inline PointContainer& point_data()
    {
        return m_points;
    };

    inline const PointContainer& point_data() const
    {
        return m_points;
    };

    inline const NodeContainer& node_data() const
    {
        return m_nodes;
    }

    inline NodeContainer& node_data()
    {
        return *m_nodes.get();
    }

    inline const IndexContainer& index_data() const
    {
        return m_indices;
    }

    inline IndexContainer& index_data()
    {
        return m_indices;
    }

    // Parameters --------------------------------------------------------------
public:
    inline int min_cell_size() const;
    inline void set_min_cell_size(int min_cell_size);

    // Internal ----------------------------------------------------------------
public:
    inline void build_rec(int node_id, int start, int end, int level);
    inline int partition(int start, int end, int dim, Scalar value);


	// Query -------------------------------------------------------------------
public :
    KdTreeKNearestPointQuery<DataPoint> k_nearest_neighbors(const VectorType& point, int k) const
    {
        return KdTreeKNearestPointQuery<DataPoint>(this, k, point);
    }

    KdTreeKNearestIndexQuery<DataPoint> k_nearest_neighbors(int index, int k) const
    {
        return KdTreeKNearestIndexQuery<DataPoint>(this, k, index);
    }

    KdTreeNearestPointQuery<DataPoint> nearest_neighbor(const VectorType& point) const
    {
        return KdTreeNearestPointQuery<DataPoint>(this, point);
    }

    KdTreeNearestIndexQuery<DataPoint> nearest_neighbor(int index) const
    {
        return KdTreeNearestIndexQuery<DataPoint>(this, index);
    }

    KdTreeRangePointQuery<DataPoint> range_neighbors(const VectorType& point, Scalar r) const
    {
        return KdTreeRangePointQuery<DataPoint>(this, r, point);
    }

    KdTreeRangeIndexQuery<DataPoint> range_neighbors(int index, Scalar r) const
    {
        return KdTreeRangeIndexQuery<DataPoint>(this, r, index);
    }

    /// Given indexs of kdTree points, support_range compute points average center and return a range  
    /// query based on this center and a sum of maximum distance of points with center and given radius.
    /// \param 
    KdTreeRangePointQuery<DataPoint> support_range_query(const vector<int>& points, Scalar radius) 
    {
        if (points.size() < 1)
        {
            std::cout << "ERROR : Points empty !" << std::endl;
            return KdTreeRangePointQuery<DataPoint>(this,0,VectorType());
        }

        VectorType center = VectorType::Zero();
        Scalar support = 0;

        for (int i = 0; i < points.size(); i++)
        {
            center = center + this->m_points[points[i]].pos();
        }

        center /= points.size();


        for (int i = 0; i < points.size(); i++)
        {
            auto t = center - this->m_points[points[i]].pos();
            Scalar temp = std::sqrt(t.squaredNorm());
            if (support < temp)
            {
                support = temp;
            }
        }

        return KdTreeRangePointQuery<DataPoint>(this,(Scalar)(support + radius) , center);
    }


    /// Given points, support_range compute points average center and return a range  
    /// query based on this center and a sum of maximum distance points-center and given radius.
    KdTreeRangePointQuery<DataPoint> support_range_query(const vector<VectorType>& points, Scalar radius) {

        if (points.size() < 1)
        {
            std::cout << "ERROR : Points empty !" << std::endl;
            return KdTreeRangePointQuery<DataPoint>(this,0,VectorType());
        }

        VectorType center = VectorType::Zero();
        Scalar support = 0;

        for (int i = 0; i < points.size(); i++)
        {
            center += points[i];
        }

        center /= points.size();


        for (int i = 0; i < points.size(); i++)
        {
            Scalar temp = std::sqrt((center - points[i]).squaredNorm());
            if (support < temp)
            {
                support = temp;
            }
        }

        return KdTreeRangePointQuery<DataPoint>(this, support + radius, center);
    }

    // Data --------------------------------------------------------------------
protected:
    PointContainer m_points;
    NodeContainer m_nodes;
    IndexContainer m_indices;

    int m_min_cell_size;
    int m_leaf_count;
};

#include "./kdTree.hpp"
} // namespace Ponca
