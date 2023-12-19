/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "enums.h"
#include "primitive.h"

#include PONCA_MULTIARCH_INCLUDE_STD(iterator)

namespace Ponca
{

#define BSKW typename BasketType::WeightFunction
#define BSKP typename BasketType::DataPoint

#ifndef PARSED_WITH_DOXYGEN
/*! \brief Namespace used for structure or classes used internally by the lib */
namespace internal
{
    template <class P, class W,
        typename Aggregate,
        template <class, class, typename> class Ext,
        template <class, class, typename> class... Exts>
    struct BasketAggregateImpl
    {
        using type = typename BasketAggregateImpl<P, W, Ext<P, W, Aggregate>, Exts...>::type;
    };

    template <class P, class W,
        typename Aggregate,
        template <class, class, typename> class Ext>
    struct BasketAggregateImpl<P, W, Aggregate, Ext>
    {
        using type = Ext<P, W, Aggregate>;
    };

    /*! \brief Internal class used to build the Basket structure */
    template <class P, class W,
        template <class, class, typename> class... Exts>
    struct BasketAggregate : BasketAggregateImpl<P, W, PrimitiveBase<P, W>, Exts...>
    {
    };

    template <typename BasketType, int Type,
        typename Aggregate,
        template <class, class, int, typename> class Ext,
        template <class, class, int, typename> class... Exts>
    struct BasketDiffAggregateImpl
    {
        using type = typename BasketDiffAggregateImpl<BasketType, Type, Ext<BSKP, BSKW, Type, Aggregate>, Exts...>::type;
    };

    template <typename BasketType, int Type,
        typename Aggregate,
        template <class, class, int, typename> class Ext>
    struct BasketDiffAggregateImpl<BasketType, Type, Aggregate, Ext>
    {
        using type = Ext<BSKP, BSKW, Type, Aggregate>;
    };

    /*! \brief Internal class used to build the BasketDiff structure */
    template <typename BasketType, int Type,
        template <class, class, int, typename> class... Exts>
    struct BasketDiffAggregate : BasketDiffAggregateImpl<BasketType, Type, BasketType, PrimitiveDer, Exts...>
    {
    };

    /*! \brief Internal class used to build the BasketAutoDiff structure */
    template <class P, class W, int Type,
        template <class, class, typename> class Ext0,
        template <class, class, typename> class... Exts>
    struct BasketAutoDiffAggregate
    {
    private:
        using Base = Ext0<P, W, PrimitiveBase<P, W>>;
        using Aggregate = typename BasketAggregateImpl<P, W, Base, Exts...>::type; // Same impl

    public:
        using type = typename Base::template DDerType<P, W, Type, PrimitiveDer<P, W, Type, Aggregate>>;
    };
}
#endif

#ifdef PONCA_CPU_ARCH
#   define WRITE_BASKET_SINGLE_HOST_FUNCTIONS                                                         \
    /*! \copydoc compute(const IteratorBegin&,const IteratorEnd&)        */                           \
    template <typename Container>                                                                     \
    PONCA_MULTIARCH inline                                                                            \
    FIT_RESULT compute(const Container& c){                                                           \
        return Self::compute(std::begin(c), std::end(c));                                             \
    }
#else
#   define WRITE_BASKET_SINGLE_HOST_FUNCTIONS
#endif

#define WRITE_BASKET_FUNCTIONS                                                                        \
    /*! \brief Convenience function for STL-like iterators               */                           \
    /*! Add neighbors stored in a container using STL-like iterators, and call finalize at the end.*/ \
    /*! The fit is evaluated multiple time if needed (see #NEED_OTHER_PASS)*/                         \
    /*! \see addNeighbor() */                                                                         \
    template <typename IteratorBegin, typename IteratorEnd>                                           \
    PONCA_MULTIARCH inline                                                                            \
    FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end){                           \
        FIT_RESULT res = UNDEFINED;                                                                   \
        do {                                                                                          \
            Self::startNewPass();                                                                     \
            for (auto it = begin; it != end; ++it){                                                   \
                Self::addNeighbor(*it);                                                               \
            }                                                                                         \
            res = Base::finalize();                                                                   \
        } while ( res == NEED_OTHER_PASS );                                                           \
        return res;                                                                                   \
    }                                                                                                 \
    /*! \brief Convenience function to iterate over a subset of samples in a PointContainer  */       \
    /*! Add neighbors stored in a PointContainer and sampled using indices stored in ids.*/           \
    /*! \tparam IndexRange STL-Like range storing indices of the neighbors */                         \
    /*! \tparam PointContainer STL-like container storing the points */                               \
    /*! \see #compute(const IteratorBegin& begin, const IteratorEnd& end)    */                       \
    template <typename IndexRange, typename PointContainer>                                           \
    PONCA_MULTIARCH inline                                                                            \
    FIT_RESULT computeWithIds(IndexRange ids, const PointContainer& points){                          \
        FIT_RESULT res = UNDEFINED;                                                                   \
        do {                                                                                          \
            Self::startNewPass();                                                                     \
            for (const auto& i : ids){                                                                \
                this->addNeighbor(points[i]);                                                         \
            }                                                                                         \
            res = this->finalize();                                                                   \
        } while ( res == NEED_OTHER_PASS );                                                           \
        return res;                                                                                   \
    }                                                                                                 \
    WRITE_BASKET_SINGLE_HOST_FUNCTIONS

    /*!
         \brief Aggregator class used to declare specialized structures with derivatives computations, using CRTP

         This is one of the central classes of the library (even if it does not perform any computation on its own).
         Checkout \ref fitting for more details, and Basket class.

         The various implementations of Ponca::Concept are mixed through specializations of the BasketDiff and Basket
         classes:
         \code
         typedef
         BasketDiff <BasketType,           // Existing Basket, to be differentiated
         DiffType,                         // Differentiation space: FitScaleDer, FitSpaceDer, or FitScaleDer|FitSpaceDer
         ComputationalDerivativesConcept1, // Implementation of ComputationalDerivativesConcept
         ComputationalDerivativesConcept2, // Implementation of ComputationalDerivativesConcept
         ... ,                             //
         > myFitDer;                       // Final structure to fit and derive a primitive over weighted samples
         \endcode

         \see Basket for the aggregation of \ref concepts_computObjectBasket "ComputationalObjectConcept"

         \tparam BasketType Existing Basket, to be differentiated
         \tparam Type Differentiation space: FitScaleDer, FitSpaceDer, or FitScaleDer|FitSpaceDer
         \tparam Ext0 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept"
         \tparam Exts Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
     */
    template <typename BasketType, int Type,
        template <class, class, int, typename> class Ext0,
        template <class, class, int, typename> class... Exts>
    class BasketDiff : public internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type {
    private:
        using Self   = BasketDiff;
    public:
    /// Base type, which aggregates all the computational objects using the CRTP
    using Base = typename internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type;
        /// Weighting function
    using WeightFunction = BSKW;
    /// Point type used for computation
    using DataPoint = BSKP;
    /// Scalar type used for computation, as defined from Basket
    using Scalar = typename DataPoint::Scalar;

    WRITE_BASKET_FUNCTIONS

    /// \copydoc Basket::addNeighbor
    PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei) {
        // compute weight
        auto wres = Base::m_w.w(_nei.pos(), _nei);
        typename Base::ScalarArray dw;

        if (wres.first > Scalar(0.)) {
            Base::addLocalNeighbor(wres.first, wres.second, _nei, dw);
            return true;
        }
        return false;
    }
};

/*!
    \brief Aggregator class used to declare specialized structures using CRTP

    This is one of the central classes of the library (even if it does not perform any computation on its own).
    Checkout \ref fitting for more details.

    The various implementations of Ponca::Concept are mixed through
    specializations of the Basket class:
    \code
        typedef
        Basket <PointImpl,           // Implementation of PointConcept
        WeightFuncImpl,              // Implementation of WeightFuncConcept
        ComputationalObjectConcept1, // Implementation of ComputationalObjectConcept
        ComputationalObjectConcept2, // Implementation of ComputationalObjectConcept
        ... ,                        //
        > myFit;                     // Final structure to fit a primitive over weighted samples
    \endcode

    \see BasketDiff for the aggregation of \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept"

    \tparam P Implements \ref ponca_concepts "PointConcept"
    \tparam W Implements \ref concepts_weighting "WeightKernelConcept"
    \tparam Ext0 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept"
    \tparam Exts Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
*/
    template <class P, class W,
        template <class, class, typename> class Ext0,
        template <class, class, typename> class... Exts>
    class Basket : public internal::BasketAggregate<P, W, Ext0, Exts...>::type
    {
    private:
        using Self   = Basket;
    public:
        /// Base type, which aggregates all the computational objects using the CRTP
        using Base = typename internal::BasketAggregate<P, W, Ext0, Exts...>::type;
        /// Scalar type used for computation, as defined from template parameter `P`
        using Scalar = typename P::Scalar;
        /// Point type used for computation
        using DataPoint = P;
        /// Weighting function
        using WeightFunction = W;

        WRITE_BASKET_FUNCTIONS;

        /// \brief Add a neighbor to perform the fit
        ///
        /// When called directly, don't forget to call PrimitiveBase::startNewPass when starting multiple passes
        /// \see compute Prefer when using a range of Points
        /// \see computeWithIds Prefer when using a range of ids
        /// \return false if param nei is not a valid neighbor (weight = 0)
        PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei) {
            // compute weight
            auto wres = Base::m_w.w(_nei.pos(), _nei);

            if (wres.first > Scalar(0.)) {
                Base::addLocalNeighbor(wres.first, wres.second, _nei);
                return true;
            }
            return false;
        }
    }; // class Basket

    template <class P, class W, int Type,
        template <class, class, typename> class Ext0,
        template <class, class, typename> class... Exts>
    class BasketAutoDiff : public internal::BasketAutoDiffAggregate<P, W, Type, Ext0, Exts...>::type
    {
    }; // class BasketAutoDiff

} //namespace Ponca

