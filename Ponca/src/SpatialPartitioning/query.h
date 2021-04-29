/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./indexSquaredDistance.h"
#include "../Common/Containers/limitedPriorityQueue.h"

#include <cmath>

namespace Ponca {

////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////

/// \brief Base class for queries storing indices
/// \ingroup spatialpartitioning
struct IndexQuery
{
    inline IndexQuery(int index = -1): m_index( index ) { }

    inline int index() const { return m_index; }
    inline void set_index(int index) { m_index = index;}

protected:
    /// Index of the queried point
    int m_index;
};


/// \brief Base class for queries storing points
/// \ingroup spatialpartitioning
template <typename DataPoint>
struct PointQuery
{
    using VectorType = typename DataPoint::VectorType;

    inline PointQuery(const VectorType& point = VectorType::Zero())
        : m_point(point) {}

    inline const VectorType& point() const { return m_point.pos(); }
    inline void set_point(const VectorType& point) { m_point = DataPoint(point); }

protected:
    /// Queried point
    DataPoint m_point;
};

/// \brief Base class for range queries
/// \ingroup spatialpartitioning
template<typename Scalar>
struct RangeQuery
{
    inline RangeQuery(Scalar radius = Scalar( 0 ))
        : m_squared_radius ( std::pow( radius, 2 ) ) {}

    inline Scalar radius() const { return std::sqrt( m_squared_radius );}
    inline Scalar squared_radius() const { return m_squared_radius; }
    inline void set_radius(Scalar radius) { m_squared_radius = std::pow( radius, 2 ); }
    inline void set_squared_radius(Scalar radius) { m_squared_radius = radius; }

protected:
	Scalar m_squared_radius { 0 };
};


/// \brief Base class for nearest queries
/// \ingroup spatialpartitioning
template<typename Scalar>
struct NearestQuery
{
    NearestQuery(){}

    int get() const{return m_nearest;}

protected:
    int m_nearest;
	Scalar m_squared_distance;
};

/// \brief Base class for knearest queries
/// \ingroup spatialpartitioning
template<typename Scalar>
struct KNearestQuery
{
    inline KNearestQuery()      : m_queue() {}
    inline KNearestQuery(int k) : m_queue(k){}

    inline size_t k() const     { return m_queue.capacity(); }
    inline void set_k(int k) { return m_queue.reserve(k); }

    inline limited_priority_queue<IndexSquaredDistance<Scalar>>& queue()
    { return m_queue; }

protected:
    limited_priority_queue<IndexSquaredDistance<Scalar>> m_queue;
};


////////////////////////////////////////////////////////////////
// KNearest Queries
////////////////////////////////////////////////////////////////

/// \brief Base class of KNearestQuery storing indices
/// \ingroup spatialpartitioning
template <typename Scalar>
struct KNearestIndexQuery : public IndexQuery,
                            public KNearestQuery<Scalar>
{
    using QueryInType  = IndexQuery;
    using QueryOutType = KNearestQuery<Scalar>;

    inline KNearestIndexQuery()
        : IndexQuery(), KNearestQuery<Scalar>() {}
    inline KNearestIndexQuery(int k)
        : IndexQuery(), KNearestQuery<Scalar>( k ) {}
    inline KNearestIndexQuery(int k, int index)
        : IndexQuery( index ), KNearestQuery<Scalar>( k ) {}
};

/// \brief Base class of KNearestQuery storing points
/// \ingroup spatialpartitioning
template <class DataPoint>
struct KNearestPointQuery : public PointQuery<DataPoint>, public KNearestQuery<typename DataPoint::Scalar>
{
    using VectorType = typename DataPoint::VectorType;
    using Scalar     = typename DataPoint::Scalar;

    using QueryInType  = PointQuery<DataPoint>;
    using QueryOutType = KNearestQuery<typename DataPoint::Scalar>;

    inline KNearestPointQuery()
        : PointQuery<DataPoint>(), KNearestQuery<Scalar>() {}
    inline KNearestPointQuery(const VectorType& point)
        : PointQuery<DataPoint>(point), KNearestQuery<Scalar>() {}
    inline KNearestPointQuery(int k)
        : PointQuery<DataPoint>(), KNearestQuery<Scalar>(k) {}
    inline KNearestPointQuery(int k, const VectorType& point)
        : PointQuery<DataPoint>(point), KNearestQuery<Scalar>(k) {}
};


////////////////////////////////////////////////////////////////
// Nearest Queries
////////////////////////////////////////////////////////////////

/// \brief Base class NearestQuery storing points
/// \ingroup spatialpartitioning
template <class DataPoint>
struct NearestPointQuery : public PointQuery<DataPoint>, public NearestQuery<typename DataPoint::Scalar>
{
	using Scalar     = typename DataPoint::Scalar;
	using VectorType = typename DataPoint::VectorType;

    using QueryInType  = PointQuery<DataPoint>;
    using QueryOutType = NearestQuery<typename DataPoint::Scalar>;

    inline NearestPointQuery()
        : PointQuery<DataPoint>(), NearestQuery<Scalar>() {}
    inline NearestPointQuery(const VectorType& point)
        : PointQuery<DataPoint>(point), NearestQuery<Scalar>() {}
};

/// \brief Base class for NearestQuery storing indices
/// \ingroup spatialpartitioning
template <typename Scalar>
struct NearestIndexQuery : public IndexQuery, public NearestQuery<Scalar>
{
    using QueryInType  = IndexQuery;
    using QueryOutType = NearestQuery<Scalar>;

    inline NearestIndexQuery() : IndexQuery(), NearestQuery<Scalar>() {}
    inline NearestIndexQuery(int index) : IndexQuery( index ), NearestQuery<Scalar>() {}
};

////////////////////////////////////////////////////////////////
// Range Queries
////////////////////////////////////////////////////////////////

/// \brief Base class RangeQuery storing points
/// \ingroup spatialpartitioning
template <class DataPoint>
struct RangePointQuery : public PointQuery<DataPoint>, public RangeQuery<typename DataPoint::Scalar>
{
	using VectorType = typename DataPoint::VectorType;
	using Scalar     = typename DataPoint::Scalar;

    using QueryInType  = PointQuery<DataPoint>;
    using QueryOutType = RangeQuery<typename DataPoint::Scalar>;

    inline RangePointQuery()
        : PointQuery<DataPoint>(), RangeQuery<Scalar>() {}
    inline RangePointQuery(Scalar radius)
        : PointQuery<DataPoint>(), RangeQuery<Scalar>(radius) {}
    inline RangePointQuery(Scalar radius, const VectorType& point)
        : PointQuery<DataPoint>(point), RangeQuery<Scalar>(radius) {}
};

/// \brief Base class RangeQuery storing indices
/// \ingroup spatialpartitioning
template<typename Scalar>
struct RangeIndexQuery : public IndexQuery, public RangeQuery<Scalar>
{
    using QueryInType  = IndexQuery;
    using QueryOutType = RangeQuery<Scalar>;

    inline RangeIndexQuery()
        : IndexQuery(), RangeQuery<Scalar>() {}
    inline RangeIndexQuery(Scalar radius)
        : IndexQuery(), RangeQuery<Scalar>( radius ) {}
    inline RangeIndexQuery(Scalar radius, int index)
        : IndexQuery( index ), RangeQuery<Scalar>( radius ) {}
};

}