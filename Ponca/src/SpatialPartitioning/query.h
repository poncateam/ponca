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
template <typename _VectorType>
struct PointQuery
{
    using VectorType = _VectorType;

    inline PointQuery(const VectorType& point = VectorType::Zero() )
        : m_point( point ) {}

    inline const VectorType& point() const { return m_point; }
    inline void set_point(const VectorType& point) { m_point = point; }

protected:
    /// Queried point
    VectorType m_point { VectorType::Zero() };
};

/// \brief Base class for range queries
/// \ingroup spatialpartitioning
struct RangeQuery
{
    inline RangeQuery(SPScalar radius = SPScalar( 0 ))
        : m_squared_radius ( std::pow( radius, 2 ) ) {}

    inline SPScalar radius() const { return std::sqrt( m_squared_radius );}
    inline SPScalar squared_radius() const { return m_squared_radius; }
    inline void set_radius(SPScalar radius) { m_squared_radius = std::pow( radius, 2 ); }
    inline void set_squared_radius(SPScalar radius) { m_squared_radius = radius; }

protected:
    SPScalar m_squared_radius { 0 };
};


/// \brief Base class for nearest queries
/// \ingroup spatialpartitioning
struct NearestQuery
{
    NearestQuery(){}

    int get() const{return m_nearest;}

protected:
    int m_nearest;
    SPScalar m_squared_distance;
};

/// \brief Base class for knearest queries
/// \ingroup spatialpartitioning
struct KNearestQuery
{
    inline KNearestQuery()      : m_queue() {}
    inline KNearestQuery(int k) : m_queue(k){}

    inline int k() const     { return m_queue.capacity(); }
    inline void set_k(int k) { return m_queue.reserve(k); }

    inline limited_priority_queue<IndexSquaredDistance>& queue()
    { return m_queue; }

protected:
    limited_priority_queue<IndexSquaredDistance> m_queue;
};


////////////////////////////////////////////////////////////////
// KNearest Queries
////////////////////////////////////////////////////////////////

/// \brief Base class of KNearestQuery storing indices
/// \ingroup spatialpartitioning
struct KNearestIndexQuery : public IndexQuery,
                            public KNearestQuery
{
    inline KNearestIndexQuery()
        : IndexQuery(), KNearestQuery() {}
    inline KNearestIndexQuery(int k)
        : IndexQuery(), KNearestQuery( k ) {}
    inline KNearestIndexQuery(int k, int index)
        : IndexQuery( index ), KNearestQuery( k ) {}
};

/// \brief Base class of KNearestQuery storing points
/// \ingroup spatialpartitioning
template <typename _VectorType>
struct KNearestPointQuery : public PointQuery<_VectorType>, public KNearestQuery
{
    using VectorType = typename PointQuery<_VectorType>::VectorType;
    inline KNearestPointQuery()
        : PointQuery<VectorType>(), KNearestQuery() {}
    inline KNearestPointQuery(const VectorType& point)
        : PointQuery<VectorType>(point), KNearestQuery() {}
    inline KNearestPointQuery(int k)
        : PointQuery<VectorType>(), KNearestQuery(k) {}
    inline KNearestPointQuery(int k, const VectorType& point)
        : PointQuery<VectorType>(point), KNearestQuery(k) {}
};


////////////////////////////////////////////////////////////////
// Nearest Queries
////////////////////////////////////////////////////////////////

/// \brief Base class NearestQuery storing points
/// \ingroup spatialpartitioning
template <typename _VectorType>
struct NearestPointQuery : public PointQuery<_VectorType>, public NearestQuery
{
    using VectorType = typename PointQuery<_VectorType>::VectorType;
    inline NearestPointQuery()
        : PointQuery<VectorType>(), NearestQuery() {}
    inline NearestPointQuery(const VectorType& point)
        : PointQuery<VectorType>(point), NearestQuery() {}
};

/// \brief Base class for NearestQuery storing indices
/// \ingroup spatialpartitioning
struct NearestIndexQuery : public IndexQuery, public NearestQuery
{
    inline NearestIndexQuery() : IndexQuery(), NearestQuery() {}
    inline NearestIndexQuery(int index) : IndexQuery( index ), NearestQuery() {}
};

////////////////////////////////////////////////////////////////
// Range Queries
////////////////////////////////////////////////////////////////

/// \brief Base class RangeQuery storing points
/// \ingroup spatialpartitioning
template <typename _VectorType>
struct  RangePointQuery : public PointQuery<_VectorType>, public RangeQuery
{
    using VectorType = typename PointQuery<_VectorType>::VectorType;

    inline RangePointQuery()
        : PointQuery<_VectorType>(), RangeQuery() {}
    inline RangePointQuery(Scalar radius)
        : PointQuery<_VectorType>(), RangeQuery( radius ) {}
    inline RangePointQuery(Scalar radius, const VectorType& point)
        : PointQuery<_VectorType>( radius ), RangeQuery( point ) {}
};

/// \brief Base class RangeQuery storing indices
/// \ingroup spatialpartitioning
struct RangeIndexQuery : public IndexQuery, public RangeQuery
{
    inline RangeIndexQuery()
        : IndexQuery(), RangeQuery() {}
    inline RangeIndexQuery(Scalar radius)
        : IndexQuery(), RangeQuery( radius ) {}
    inline RangeIndexQuery(Scalar radius, int index)
        : IndexQuery( index ), RangeQuery( radius ) {}
};

}