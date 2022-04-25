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


/// \internal
/// \brief Macro generating code of the the Query base classes inhering QueryInputIsIndex
/// \note For internal use only
#define DECLARE_INDEX_QUERY_CLASS(OUT_TYPE) \
/*! \brief Base Query class combining QueryInputIsIndex and QueryOutputIs##OUT_TYPE##. */    \
template <typename Scalar>                                 \
struct  OUT_TYPE##IndexQuery : Query<QueryInputIsIndex, QueryOutputIs##OUT_TYPE<Scalar>> \
{ \
    using Base = Query<QueryInputIsIndex, QueryOutputIs##OUT_TYPE<Scalar>>; \
    using Base::Base; \
};


/// \internal
/// \brief Macro generating code of the the Query base classes inhering QueryInputIsPosition
/// \note For internal use only
#define DECLARE_POINT_QUERY_CLASS(OUT_TYPE) \
/*! \brief Base Query class combining QueryInputIsPosition and QueryOutputIs##OUT_TYPE##. */    \
template <typename DataPoint>                                \
struct  OUT_TYPE##PointQuery : Query<QueryInputIsPosition<DataPoint>, \
                                     QueryOutputIs##OUT_TYPE<typename DataPoint::Scalar>> \
{ \
    using Base = Query<QueryInputIsPosition<DataPoint>, QueryOutputIs##OUT_TYPE<typename DataPoint::Scalar>>; \
    using Base::Base; \
};

/// \addtogroup spatialpartitioning
/// @{

////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////

/// \brief Base class for queries input type
    struct QueryInputBase {
    };

/// \brief Base class for queries output type
    struct QueryOutputBase {
        struct DummyOutputParameter {
        };
    };

/// \brief Base class for typed queries input type
    template<typename InputType_>
    struct QueryInput : public QueryInputBase {
        using InputType = InputType_;

        inline QueryInput(InputType input) : m_input(input) {}

        inline const InputType &input() const { return m_input; }
    protected:
        /// \brief Edit the input (queries need to be restarted for good behavior)
        ///Usefull to avoid query reallocation between different requests
        inline void editInput(const InputType& input) { m_input = input; }
    
    private:
        /// Index of the queried point
        const InputType m_input;
    };


/// \brief Base class for queries storing points
    struct QueryInputIsIndex : public QueryInput<int> {
        using Base = QueryInput<int>;
        using InputType = typename Base::InputType;

        inline QueryInputIsIndex(const InputType &point = -1)
                : Base(point) {}
    };

/// \brief Base class for queries storing points
    template<typename DataPoint>
    struct QueryInputIsPosition : public QueryInput<typename DataPoint::VectorType> {
        using Base = QueryInput<typename DataPoint::VectorType>;
        using InputType = typename Base::InputType;

        inline QueryInputIsPosition(const InputType &point = InputType::Zero())
                : Base(point) {}
    };

/// \brief Base class for range queries
    template<typename Scalar>
    struct QueryOutputIsRange : public QueryOutputBase {
        using OutputParameter = Scalar;

        inline QueryOutputIsRange(OutputParameter radius = OutputParameter(0))
                : m_squared_radius(std::pow(radius, 2)) {}

        inline Scalar radius() const { return std::sqrt(m_squared_radius); }

        inline Scalar squared_radius() const { return m_squared_radius; }

        inline void set_radius(Scalar radius) { m_squared_radius = std::pow(radius, 2); }

        inline void set_squared_radius(Scalar radius) { m_squared_radius = radius; }

    protected:
        /// \brief Reset Query for a new search
        inline void reset() { }
        Scalar m_squared_radius{0};
    };

/// \brief Base class for nearest queries
    template<typename Scalar>
    struct QueryOutputIsNearest : public QueryOutputBase {
        using OutputParameter = typename QueryOutputBase::DummyOutputParameter;

        QueryOutputIsNearest() {}

        int get() const { return m_nearest; }

    protected:
        /// \brief Reset Query for a new search
        void reset() {
            m_nearest = -1;
            m_squared_distance = std::numeric_limits<Scalar>::max();
        }

        int m_nearest {-1};
        Scalar m_squared_distance {std::numeric_limits<Scalar>::max()};
    };

/// \brief Base class for knearest queries
    template<typename Scalar>
    struct QueryOutputIsKNearest : public QueryOutputBase {
        using OutputParameter = int;

        inline QueryOutputIsKNearest(OutputParameter k = 0) : m_queue(k) {}

        inline limited_priority_queue<IndexSquaredDistance<Scalar>> &queue() { return m_queue; }

    protected:
        /// \brief Reset Query for a new search
        void reset() {
            m_queue.clear();
            m_queue.push({-1,std::numeric_limits<Scalar>::max()});
        }
        limited_priority_queue<IndexSquaredDistance<Scalar>> m_queue;
    };


    template<typename Input_, typename Output_>
    struct Query : public Input_, public Output_ {
        using QueryInType = Input_;
        using QueryOutType = Output_;

        static_assert(std::is_base_of<QueryInputBase, QueryInType>::value,
                      "QueryInType must inherit Ponca::QueryInputBase");
        static_assert(std::is_base_of<QueryOutputBase, QueryOutType>::value,
                      "QueryInType must inherit Ponca::QueryInputBase");

        inline Query(const typename QueryInType::InputType &in)
                : QueryInType(in), QueryOutType() {}

        inline Query(const typename QueryOutType::OutputParameter &outParam,
                     const typename QueryInType::InputType &in)
                : QueryOutType(outParam), QueryInType(in) {}
    };

DECLARE_INDEX_QUERY_CLASS(KNearest) //KNearestIndexQuery
DECLARE_INDEX_QUERY_CLASS(Nearest)  //NearestIndexQuery
DECLARE_INDEX_QUERY_CLASS(Range)    //RangeIndexQuery
DECLARE_POINT_QUERY_CLASS(KNearest) //KNearestPointQuery
DECLARE_POINT_QUERY_CLASS(Nearest)  //NearestPointQuery
DECLARE_POINT_QUERY_CLASS(Range)    //RangePointQuery

/// @}

#undef DECLARE_INDEX_QUERY_CLASS
#undef DECLARE_POINT_QUERY_CLASS
}