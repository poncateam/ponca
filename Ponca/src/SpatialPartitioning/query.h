/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./indexSquaredDistance.h"
#include "../Common/Containers/limitedPriorityQueue.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

namespace Ponca {


/// \internal
/// \brief Macro generating code of the the Query base classes inhering QueryInputIsIndex
/// \note For internal use only
#define DECLARE_INDEX_QUERY_CLASS(OUT_TYPE) \
/*! \brief Base Query class combining QueryInputIsIndex and QueryOutputIs##OUT_TYPE##. */    \
template <typename Index, typename Scalar> \
struct  OUT_TYPE##IndexQuery : Query<QueryInputIsIndex<Index>, QueryOutputIs##OUT_TYPE<Index, Scalar>> \
{ \
    using Base = Query<QueryInputIsIndex<Index>, QueryOutputIs##OUT_TYPE<Index, Scalar>>; \
    using Base::Base; \
};


/// \internal
/// \brief Macro generating code of the the Query base classes inhering QueryInputIsPosition
/// \note For internal use only
#define DECLARE_POINT_QUERY_CLASS(OUT_TYPE) \
/*! \brief Base Query class combining QueryInputIsPosition and QueryOutputIs##OUT_TYPE##. */    \
template <typename Index, typename DataPoint> \
struct  OUT_TYPE##PointQuery : Query<QueryInputIsPosition<DataPoint>, \
                                     QueryOutputIs##OUT_TYPE<Index, typename DataPoint::Scalar>> \
{ \
    using Base = Query<QueryInputIsPosition<DataPoint>, QueryOutputIs##OUT_TYPE<Index, typename DataPoint::Scalar>>; \
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
        /// \brief Edit the input 
        /// Need to be used carefully. Modifying a query input while iterating on the query will result in undefined behavior.
        /// Simplest way to avoid this is to restart the iteration on the query. 
        /// Usefull to avoid query reallocation between different requests
        inline void setInput(const InputType& input) { m_input = input; }
    
    private:
        /// Index of the queried point
        InputType m_input;
    };


/// \brief Base class for queries storing points
    template <typename Index>
    struct QueryInputIsIndex : public QueryInput<Index> {
        using Base = QueryInput<Index>;
        using InputType = typename Base::InputType;

        inline QueryInputIsIndex(const InputType &point = -1)
                : Base(point) {}

        inline void operator()(const InputType &point = InputType::Zero()){
            Base::setInput(point);
        }
    protected:
        /// Functor used to check if a given Idx must be skipped
        template <typename IndexType>
        inline bool skipIndexFunctor(IndexType idx) const {return Base::input() == idx;};
        /// Generic method to access input position. Container is expected to hold kdtree positions
        template <typename Container>
        inline auto getInputPosition(const Container &c) -> const typename Container::value_type::VectorType
        { return c[Base::input()].pos(); }
    };

/// \brief Base class for queries storing points
    template<typename DataPoint>
    struct QueryInputIsPosition : public QueryInput<typename DataPoint::VectorType> {
        using Base = QueryInput<typename DataPoint::VectorType>;
        using InputType = typename Base::InputType;

        inline QueryInputIsPosition(const InputType &point = InputType::Zero())
                : Base(point) {}

        inline void operator()(const InputType &point = InputType::Zero()){
            Base::setInput( point );
        }
    protected:
        /// Functor used to check if a given Idx must be skipped
        template <typename IndexType>
        inline bool skipIndexFunctor(IndexType idx) const {return false;};
        /// Generic method to access input position. Container is expected to hold kdtree positions
        template <typename Container>
        inline auto getInputPosition(const Container &) -> const typename Container::value_type::VectorType
        { return Base::input(); }
    };

/// \brief Base class for range queries
    template<typename Index, typename Scalar>
    struct QueryOutputIsRange : public QueryOutputBase {
        using OutputParameter = Scalar;

        inline QueryOutputIsRange(OutputParameter radius = OutputParameter(0))
                : m_squared_radius(PONCA_MULTIARCH_CU_STD_NAMESPACE(pow)(radius, OutputParameter(2))) {}

        inline void operator() (OutputParameter radius){
            setRadius( radius );
        }

        inline Scalar radius() const {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            return sqrt(m_squared_radius);
        }

        inline Scalar squaredRadius() const { return m_squared_radius; }

        inline void setRadius(Scalar radius) {
            setSquaredRadius (radius*radius);
        }

        inline void setSquaredRadius(Scalar radius) { m_squared_radius = radius; }

    protected:
        /// \brief Reset Query for a new search
        inline void reset() { }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        inline Scalar descentDistanceThreshold() const { return m_squared_radius; }
        Scalar m_squared_radius{0};
    };

/// \brief Base class for nearest queries
    template<typename Index, typename Scalar>
    struct QueryOutputIsNearest : public QueryOutputBase {
        using OutputParameter = typename QueryOutputBase::DummyOutputParameter;

        QueryOutputIsNearest() {}

        inline void operator() (){ }

        Index get() const { return m_nearest; }

    protected:
        /// \brief Reset Query for a new search
        void reset() {
            PONCA_MULTIARCH_STD_MATH(numeric_limits);
            m_nearest = -1;
            m_squared_distance = numeric_limits<Scalar>::max();
        }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        inline Scalar descentDistanceThreshold() const { return m_squared_distance; }

        Index m_nearest {-1};
        Scalar m_squared_distance {PONCA_MULTIARCH_CU_STD_NAMESPACE(numeric_limits)<Scalar>::max()};
    };

/// \brief Base class for knearest queries
    template<typename Index, typename Scalar>
    struct QueryOutputIsKNearest : public QueryOutputBase {
        using OutputParameter = Index;

        inline QueryOutputIsKNearest(OutputParameter k = 0) : m_queue(k) {}

        inline void operator() (OutputParameter k) { m_queue = limited_priority_queue<IndexSquaredDistance<Index, Scalar>>(k); }

        inline limited_priority_queue<IndexSquaredDistance<Index, Scalar>> &queue() { return m_queue; }

    protected:
        /// \brief Reset Query for a new search
        void reset() {
            PONCA_MULTIARCH_STD_MATH(numeric_limits);
            m_queue.clear();
            m_queue.push({-1, numeric_limits<Scalar>::max()});
        }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        inline Scalar descentDistanceThreshold() const { return m_queue.bottom().squared_distance; }
        limited_priority_queue<IndexSquaredDistance<Index, Scalar>> m_queue;
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

        template<typename Base, typename... outputType>
        inline Base& operator()(const typename QueryInType::InputType &in, outputType&&... out){
            QueryInType:: operator()(in);
            QueryOutType::operator()(std::forward<outputType>(out)...);
            return *((Base*)(this));
        }

        template<typename Base>
        inline Base& operator()(const typename QueryInType::InputType &in){
            QueryInType:: operator()(in);
            return *((Base*)(this));
        }
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