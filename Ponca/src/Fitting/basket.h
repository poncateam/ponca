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
}
#endif

    /*!
         \brief ComputeObject is a virtual object that represents an algorithm which can be used with the compute functions.
         The compute(begin, end) and computeWithIds(ids, points) methods must be implemented by the inheriting class.
         \note The compute(container) that is defined in this structure can be reused in the inheriting class by adding
         "using ComputeObject<Self>::compute;" to make it accessible
    */
    template <typename Derived>
    struct ComputeObject {
    protected:
        /// \brief Retrieve the top layer object
        /// Returns a reference to the derived class so that we can use its overwritten methods
        Derived& derived() { return static_cast<Derived&>(*this); }
    public:

#ifdef PONCA_CPU_ARCH
        /*! \brief Convenience function for STL-like container                                                        */
        /*! Uses the compute(IteratorBegin, IteratorEnd) function                                                     */
        /*! \note This method is only accessible when using a CPU architecture (PONCA_CPU_ARCH = true)                */
        /*! \tparam Container And STL-Like container                                                                  */
        /*! \see #compute(const IteratorBegin& begin, const IteratorEnd& end)                                         */
        template <typename Container>
        FIT_RESULT compute(const Container& c) {
            return derived().compute(std::begin(c), std::end(c));
        }
#endif

        /*! \brief Convenience function for STL-like iterators
            \tparam IteratorBegin The beginning of the iterator (std::begin(iterator)
            \tparam IteratorEnd   The end of the iterator (std::end(iterator)
        */
        template <typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end)
        {
            std::cerr << "ERROR" << std::endl;
            return UNDEFINED;
        };

        /*! \brief Convenience function to iterate over a subset of samples in a PointContainer
            \tparam IndexRange STL-Like range storing indices of the neighbors
            \tparam PointContainer STL-like container storing the points
            \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
        */
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(IndexRange /*ids*/, const PointContainer& /*points*/) {
            std::cerr << "ERROR" << std::endl;
            return UNDEFINED;
        };
    }; // struct ComputeObject

    template<typename Derived, typename Base>
    struct BasketComputeObject : public ComputeObject<Derived> {
    protected:
        using ComputeObject<Derived>::derived;
        Base& base() { return static_cast<Base&>(static_cast<Derived&>(*this)); }
    public:
        using ComputeObject<Derived>::compute; // Makes the default compute(container) accessible when using a CPU architecture

        /*!
         * \brief Convenience function for STL-like iterators
         * Add neighbors stored in a container using STL-like iterators, and call finalize at the end.
         * The fit is evaluated multiple time if needed (see #NEED_OTHER_PASS)
         * \see addNeighbor()
         */
        template <typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end){
            base().init();
            FIT_RESULT res = UNDEFINED;

            do {
                derived().startNewPass();
                for (auto it = begin; it != end; ++it){
                    derived().addNeighbor(*it);
                }
                res = base().finalize();
            } while ( res == NEED_OTHER_PASS );

            return res;
        }

        /*!
         * \brief Convenience function to iterate over a subset of samples in a PointContainer
         * Add neighbors stored in a PointContainer and sampled using indices stored in ids.
         * \tparam IndexRange STL-Like range storing indices of the neighbors
         * \tparam PointContainer STL-like container storing the points
         * \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
         */
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(IndexRange ids, const PointContainer& points){
            base().init();
            FIT_RESULT res = UNDEFINED;

            do {
                derived().startNewPass();
                for (const auto& i : ids){
                    derived().addNeighbor(points[i]);
                }
                res = base().finalize();
            } while ( res == NEED_OTHER_PASS );

            return res;
        }
    };
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
    class BasketDiff : public BasketComputeObject<BasketDiff<BasketType, Type, Ext0, Exts...>, typename internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type>,
                       public internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type {
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

    using BasketComputeObject<Self, Base>::compute;
    using BasketComputeObject<Self, Base>::computeWithIds;

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
    class Basket : public BasketComputeObject<Basket<P, W, Ext0, Exts...>, typename internal::BasketAggregate<P, W, Ext0, Exts...>::type>,
                   public internal::BasketAggregate<P, W, Ext0, Exts...>::type
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

        using BasketComputeObject<Self, Base>::compute;
        using BasketComputeObject<Self, Base>::computeWithIds;

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

} //namespace Ponca

