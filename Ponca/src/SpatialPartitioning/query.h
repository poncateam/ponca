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
/// \brief Macro generating code of the Query base classes inhering QueryInputIsIndex.
/// \note For internal use only
#define DECLARE_INDEX_QUERY_CLASS(OUT_TYPE)                                                                            \
/*! \brief Base Query class combining QueryInputIsIndex and QueryOutputIs##OUT_TYPE##.                              */ \
/*! `IndexQuery` objects acts as a `Range` that can be iterated over.                                               */ \
/*! They are used as the return type for the index searches                                                         */ \
/*! and they allow easy access to the result outputs.                                                               */ \
/*! This specialization of the `IndexQuery` concept is used to iterate over the neighbors of a given point index,   */ \
/*! using a ##OUT_TYPE## Index Query request.                                                                       */ \
template <typename Index, typename Scalar>                                                                             \
struct  OUT_TYPE##IndexQuery : Query<QueryInputIsIndex<Index>, QueryOutputIs##OUT_TYPE<Index, Scalar>>                 \
{                                                                                                                      \
    /*! \brief Alias to the inherited base Query  */                                                                   \
    using Base = Query<QueryInputIsIndex<Index>, QueryOutputIs##OUT_TYPE<Index, Scalar>>;                              \
    /*! \brief Inherited default constructor  */                                                                       \
    using Base::Base;                                                                                                  \
};


/// \internal
/// \brief Macro generating code of the Query base classes inhering QueryInputIsPosition.
/// \note For internal use only
#define DECLARE_POINT_QUERY_CLASS(OUT_TYPE)                                                                            \
/*! \brief Base Query class combining QueryInputIsPosition and QueryOutputIs##OUT_TYPE##.                           */ \
/*! `PointQuery` objects acts as a `Range` that can be iterated over.                                               */ \
/*! They are used as the return type for the index searches                                                         */ \
/*! and they allow easy access to the result outputs.                                                               */ \
/*! This specialization of the `PointQuery` concept is used to iterate over the neighbors of a given point position,*/ \
/*! using a ##OUT_TYPE## Point Query request.                                                                       */ \
template <typename Index, typename DataPoint>                                                                          \
struct  OUT_TYPE##PointQuery : Query<QueryInputIsPosition<DataPoint>,                                                  \
                                     QueryOutputIs##OUT_TYPE<Index, typename DataPoint::Scalar>>                       \
{                                                                                                                      \
    /*! \brief Alias to the inherited base Query  */                                                                   \
    using Base = Query<QueryInputIsPosition<DataPoint>, QueryOutputIs##OUT_TYPE<Index, typename DataPoint::Scalar>>;   \
    /*! \brief Inherited default constructor  */                                                                       \
    using Base::Base;                                                                                                  \
};

/// \addtogroup spatialpartitioning
/// @{

////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////

/// \brief Base class for queries input type
    struct QueryInputBase {
    };

    /*!
     * \brief Base class for queries output types.
     *
     * `QueryOutput` objects are the return types of the searches of the accelerating structures
     *  defined in the SpatialPartitioning module.
     *
     * They are used for easy access to the result outputs.
     * (e.g. to iterate over the neighbors of the evaluated point using a `rangeNeighborsQuery` request).
     */
    struct QueryOutputBase {
        struct DummyOutputParameter {
        };
    };

    /*!
     * \brief Base class for Query input type.
     *
     * Stores internally a value that is related to the point about to be evaluated
     * (e.g. the position of the point, or its index)
     *
     * \warning This class has to be specialized for a specific input type, and can't be used as is.
     * \see QueryInputIsIndex, QueryInputIsPosition
     */
    template<typename InputType_>
    struct QueryInput : public QueryInputBase {
        /// \brief Alias to the templated input type
        using InputType = InputType_;

        /// \brief Default constructor that initialize the input field
        PONCA_MULTIARCH inline QueryInput(InputType input) : m_input(input) {}

        /// \brief Read access to the input
        PONCA_MULTIARCH inline const InputType &input() const { return m_input; }
    protected:
        /*!
         * \brief Edit the input
         *
         * \warning Need to be used carefully.
         * Modifying a query input while iterating on the query will result in undefined behavior.
         * Simplest way to avoid this is to restart the iteration on the query.
         * Usefull to avoid query reallocation between different requests.
         */
        PONCA_MULTIARCH inline void setInput(const InputType& input) { m_input = input; }
    
    private:
        /// Index of the queried point
        InputType m_input;
    };


    /*!
     * \brief Extension of `QueryInput` that handles an **index** based search, in a partitioning structure.
     *
     * Stores internally the index of the evaluated point.
     * \see QueryInput
     */
    template <typename Index>
    struct QueryInputIsIndex : public QueryInput<Index> {
        /// \brief Alias to the inherited type
        using Base = QueryInput<Index>;
        /// \brief Alias to the templated input type (Index)
        using InputType = typename Base::InputType;

        /// \brief Default constructor that initialize the input index
        PONCA_MULTIARCH inline QueryInputIsIndex(const InputType &point = -1)
                : Base(point) {}

        /// \brief Access operator that resets the input index
        PONCA_MULTIARCH inline void operator()(const InputType &point = InputType::Zero()){
            Base::setInput(point);
        }
    protected:
        /// Functor used to check if a given Idx must be skipped
        template <typename IndexType>
        PONCA_MULTIARCH inline bool skipIndexFunctor(IndexType idx) const {return Base::input() == idx;};
        /// Generic method to access input position. The VectorType needs to be passed as a template argument.
        template <typename VectorType, typename Container>
        PONCA_MULTIARCH inline const VectorType getInputPosition(const Container &c)
        { return c[Base::input()].pos(); }
    };

    /*!
     * \brief Extension of `QueryInput` that handles a **position** based search, in a partitioning structure.
     *
     * Stores internally the position of the evaluated point.
     * \see QueryInput
     */
    template<typename DataPoint>
    struct QueryInputIsPosition : public QueryInput<typename DataPoint::VectorType> {
        /// \brief Alias to the inherited type
        using Base = QueryInput<typename DataPoint::VectorType>;
        /// \brief Alias to the templated input type (Index)
        using InputType = typename Base::InputType;

        /// \brief Default constructor that initialize the input position
        PONCA_MULTIARCH inline QueryInputIsPosition(const InputType &point = InputType::Zero())
                : Base(point) {}

        /// \brief Access operator that resets the input point
        PONCA_MULTIARCH inline void operator()(const InputType &point = InputType::Zero()){
            Base::setInput( point );
        }
    protected:
        /// Functor used to check if a given Idx must be skipped
        template <typename IndexType>
        PONCA_MULTIARCH inline bool skipIndexFunctor(IndexType idx) const {return false;};
        /// Generic method to access input position. The VectorType needs to be passed as a template argument.
        template <typename VectorType=typename DataPoint::VectorType, typename Container>
        PONCA_MULTIARCH inline const VectorType getInputPosition(const Container &)
        { return Base::input(); }
    };

    /*!
     * \brief Class to construct the range query output.
     *
     * Stores internally the radius value of a range request.
     * \see QueryOutputBase
     */
    template<typename Index, typename Scalar>
    struct QueryOutputIsRange : public QueryOutputBase {
        /// \brief Alias to Output type
        using OutputParameter = Scalar;

        /// \brief Default constructor that initialize the output parameter value
        PONCA_MULTIARCH inline QueryOutputIsRange(OutputParameter radius = OutputParameter(0))
                : m_squared_radius(PONCA_MULTIARCH_STD_MATH_NAMESPACE(pow)(radius, OutputParameter(2))) {}

        /// \brief Access operator that resets the output parameter
        PONCA_MULTIARCH inline void operator() (OutputParameter radius){
            setRadius( radius );
        }

        /*!
         * \brief Generic method to access the radius.
         *
         * \warning This getter method is more expensive to process than `squaredRadius`,
         * because it computes the square root of the squared radius each time it is called.
         * This is done this way to avoid having to store internally two values related to the radius.
         *
         * `squaredRadius` is overall better for distance comparison.
         * \see squaredRadius
         */
        PONCA_MULTIARCH inline Scalar radius() const {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            return sqrt(m_squared_radius);
        }

        /// \brief Generic method to access the radius squared.
        PONCA_MULTIARCH inline Scalar squaredRadius() const { return m_squared_radius; }

        /*!
         * \brief Set the radius distance of the query
         *
         * \note Store internally the squared radius for faster distance comparison
         */
        PONCA_MULTIARCH inline void setRadius(Scalar radius) {
            setSquaredRadius (radius*radius);
        }

        /// \brief Set the squared radius distance of the query
        PONCA_MULTIARCH inline void setSquaredRadius(Scalar radius) { m_squared_radius = radius; }

    protected:
        /// \brief Reset Query for a new search
        PONCA_MULTIARCH inline void reset() { }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        PONCA_MULTIARCH inline Scalar descentDistanceThreshold() const { return m_squared_radius; }
        /// \brief Radius used for the search
        Scalar m_squared_radius{0};
    };

    /*!
     * \brief Class to construct the nearest query output.
     *
     * Stores internally the nearest neighbor and the Distance threshold (for tree descent).
     * \see QueryOutputBase
     */
    template<typename Index, typename Scalar>
    struct QueryOutputIsNearest : public QueryOutputBase {
        /// \brief Alias to Output type
        using OutputParameter = typename QueryOutputBase::DummyOutputParameter;

        /// \brief Default constructor
        PONCA_MULTIARCH QueryOutputIsNearest() {}

        /// \brief Access operator
        PONCA_MULTIARCH inline void operator() (){ }

        /// \brief Get the closest points
        PONCA_MULTIARCH Index get() const { return m_nearest; }

    protected:
        /// \brief Reset Query for a new search
        PONCA_MULTIARCH void reset() {
            m_nearest = -1;
            m_squared_distance = PONCA_MULTIARCH_CU_STD_NAMESPACE(numeric_limits)<Scalar>::max();
        }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        PONCA_MULTIARCH inline Scalar descentDistanceThreshold() const { return m_squared_distance; }

        /// \brief Index of the nearest neighbor
        Index m_nearest {-1};
        /// \brief Distance to the nearest neighbor
        Scalar m_squared_distance {PONCA_MULTIARCH_CU_STD_NAMESPACE(numeric_limits)<Scalar>::max()};
    };

    /*! \brief Class to construct the knearest queries
     *
     *  Stores internally the neighbors collection of the knn request and the Distance threshold (for tree descent).
     *  \see QueryOutputBase
     */
    template<typename Index, typename Scalar>
    struct QueryOutputIsKNearest : public QueryOutputBase {
        /// \brief Alias to Output type
        using OutputParameter = Index;

        /// \brief Default constructor that initialize the output parameter value
        PONCA_MULTIARCH inline QueryOutputIsKNearest(OutputParameter k = 0) : m_queue(k) {}

        /// \brief Access operator that resets the output parameter
        PONCA_MULTIARCH inline void operator() (OutputParameter k) { m_queue = limited_priority_queue<IndexSquaredDistance<Index, Scalar>>(k); }

        /// \brief Access to the priority queue storing the neighbors
        PONCA_MULTIARCH inline limited_priority_queue<IndexSquaredDistance<Index, Scalar>> &queue() { return m_queue; }

    protected:
        /// \brief Reset Query for a new search
        PONCA_MULTIARCH void reset() {
            m_queue.clear();
            m_queue.push({-1, PONCA_MULTIARCH_CU_STD_NAMESPACE(numeric_limits)<Scalar>::max()});
        }
        /// \brief Distance threshold used during tree descent to select nodes to explore
        PONCA_MULTIARCH inline Scalar descentDistanceThreshold() const { return m_queue.bottom().squared_distance; }
        /// \brief Queue storing the neighbors
        limited_priority_queue<IndexSquaredDistance<Index, Scalar>> m_queue;
    };

    /*!
     * \brief Composes the Query object depending on an input type and output type.
     *
     * \tparam Input_ The query input type corresponds to the value used for the search.
     * \tparam Output_ The query output type corresponds to the results of the search (it can be iterated over).
     */
    template<typename Input_, typename Output_>
    struct Query : public Input_, public Output_ {
        /// \brief Alias to the input type
        using QueryInType = Input_;
        /// \brief Alias to the output type
        using QueryOutType = Output_;

        static_assert(std::is_base_of<QueryInputBase, QueryInType>::value,
                      "QueryInType must inherit Ponca::QueryInputBase");
        static_assert(std::is_base_of<QueryOutputBase, QueryOutType>::value,
                      "QueryInType must inherit Ponca::QueryInputBase");

        /// \brief Default constructor that initialize the input parameter only
        PONCA_MULTIARCH inline Query(const typename QueryInType::InputType &in)
                : QueryInType(in), QueryOutType() {}

        /// \brief Default constructor that initialize the input and output parameters
        PONCA_MULTIARCH inline Query(const typename QueryOutType::OutputParameter &outParam,
                     const typename QueryInType::InputType &in)
                : QueryOutType(outParam), QueryInType(in) {}

        /// \brief Access operator that resets the input and output parameters
        template<typename Base, typename... outputType>
        PONCA_MULTIARCH inline Base& operator()(const typename QueryInType::InputType &in, outputType&&... out){
            QueryInType:: operator()(in);
            QueryOutType::operator()(std::forward<outputType>(out)...);
            return *((Base*)(this));
        }

        /// \brief Access operator that resets the input parameter only
        template<typename Base>
        PONCA_MULTIARCH inline Base& operator()(const typename QueryInType::InputType &in){
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