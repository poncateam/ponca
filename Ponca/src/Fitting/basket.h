/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "enums.h"
#include "primitive.h"
#include "compute.h"

#include PONCA_MULTIARCH_INCLUDE_STD(iterator)

namespace Ponca
{

#define BSKW typename BasketType::WeightFunction
#define BSKP typename BasketType::DataPoint

#ifndef PARSED_WITH_DOXYGEN
/*! \brief Namespace used for structure or classes used internally by the lib */
namespace internal
{
    /*!
     * \brief This class unrolls the extension (from left to right) from the CRTP variadic list
     *
     * \tparam P Implements \ref ponca_concepts "PointConcept"
     * \tparam W Implements \ref concepts_weighting "WeightKernelConcept"
     * \tparam Aggregate The base CRTP class
     * \tparam Ext First (left-side) extension in the CTRP variadic list
     * \tparam Exts Remaining elements (excluding Ext) of the CRTP variadic list
     */
    template <class P, class W,
        typename Aggregate,
        template <class, class, typename> class Ext,
        template <class, class, typename> class... Exts>
    struct BasketAggregateImpl
    {
        using type = typename BasketAggregateImpl<P, W, Ext<P, W, Aggregate>, Exts...>::type;
    };

    /*!
     * \brief Specialized version of BasketAggregateImpl when the variadic list is empty
     *
     * \tparam P Implements \ref ponca_concepts "PointConcept"
     * \tparam W Implements \ref concepts_weighting "WeightKernelConcept"
     * \tparam Aggregate The base CRTP class
     * \tparam Ext Unique (or last) extension of the CTRP variadic list
     */
    template <class P, class W,
        typename Aggregate,
        template <class, class, typename> class Ext>
    struct BasketAggregateImpl<P, W, Aggregate, Ext>
    {
        using type = Ext<P, W, Aggregate>;
    };

    /*!
     * \brief Internal class used to build the Basket structure
     * Uses BasketAggregateImpl to unroll the CRTP variadic list
     *
     * \tparam P Implements \ref ponca_concepts "PointConcept"
     * \tparam W Implements \ref concepts_weighting "WeightKernelConcept"
     * \tparam Exts CRTP variadic list
     */
    template <class P, class W,
        template <class, class, typename> class... Exts>
    struct BasketAggregate : BasketAggregateImpl<P, W, PrimitiveBase<P, W>, Exts...>
    {
    };

    /*!
     * \brief This class unrolls the extension (from left to right) from the CRTP variadic list
     * \see BasketDiffAggregateImpl
     *
     * \tparam Type Differentiation type
     * \tparam BasketType Input Basket Type
     * \tparam Ext First (left-side) extension in the CTRP variadic list
     * \tparam Exts Remaining elements (excluding Ext) of the CRTP variadic list
     */
    template <int Type,
        typename BasketType,
        template <class, class, int, typename> class Ext,
        template <class, class, int, typename> class... Exts>
    struct BasketDiffAggregateImpl
    {
        using type = typename BasketDiffAggregateImpl<Type, Ext<BSKP, BSKW, Type, BasketType>, Exts...>::type;
    };

    /*!
     * \brief Specialized version of BasketDiffAggregateImpl when the variadic list is empty
     *
     * \tparam Type Differentiation type
     * \tparam BasketType Input Basket Type
     * \tparam Ext Unique (or last) extension of the CTRP variadic list
     */
    template <int Type,
        typename BasketType,
        template <class, class, int, typename> class Ext>
    struct BasketDiffAggregateImpl<Type, BasketType, Ext>
    {
        using type = Ext<BSKP, BSKW, Type, BasketType>;
    };

    /*!
     * \brief Internal class used to build the BasketDiff structure
     * Uses BasketDiffAggregateImpl to unroll the CRTP variadic list
     *
     * \tparam BasketType BasketType Input Basket Type
     * \tparam Type Differentiation type
     * \tparam Exts CRTP variadic list
     */
    template <typename BasketType, int Type,
        template <class, class, int, typename> class... Exts>
    struct BasketDiffAggregate : BasketDiffAggregateImpl<Type, BasketType, PrimitiveDer, Exts...>
    {
    };
}
#endif

    /*!
     * Base ComputeObject for the Basket classes
     *
     * Implements the compute methods for fitting: #compute, #computeWithIds, ...
     * Checkout \ref fitting for more details
     *
     * The various implementations of Ponca::Concept are mixed through specializations of the BasketDiff and Basket
     *   classes:
     *   \code
     *   typedef
     *   BasketDiff <BasketType,           // Existing Basket, to be differentiated
     *   DiffType,                         // Differentiation space: FitScaleDer, FitSpaceDer, or FitScaleDer|FitSpaceDer
     *   ComputationalDerivativesConcept1, // Implementation of ComputationalDerivativesConcept
     *   ComputationalDerivativesConcept2, // Implementation of ComputationalDerivativesConcept
     *   ... ,                             //
     *   > myFitDer;                       // Final structure to fit and derive a primitive over weighted samples
     *   \endcode
     *
     * @tparam Derived Derived class that provides the addNeighbor method (either Basket or BasketDiff)
     * @tparam Base Base class that provides, through the CRTP the init, startNewPass, addNeighbor and finalize methods
     */
    template<typename _Derived, typename _Base>
    struct BasketComputeObject : public ComputeObject<_Derived>, public virtual _Base {
        using Base    = _Base;    /// <\brief Alias to the Base type
        using Derived = _Derived; /// \brief Alias to the Derived type
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

#define WRITE_COMPUTE_FUNCTIONS \
    using BasketComputeObject<Self, Base>::compute; \
    using BasketComputeObject<Self, Base>::computeWithIds;

    /*!
         \brief Aggregator class used to declare specialized structures with derivatives computations, using CRTP
         \copydoc BasketComputeObject
         \tparam BasketType Existing Basket, to be differentiated
         \tparam Type Differentiation space: FitScaleDer, FitSpaceDer, or FitScaleDer|FitSpaceDer
         \tparam Ext0 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept"
         \tparam Exts Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
     */
    template <typename BasketType, int Type,
        template <class, class, int, typename> class Ext0,
        template <class, class, int, typename> class... Exts>
    class BasketDiff : public BasketComputeObject<BasketDiff<BasketType, Type, Ext0, Exts...>,
                                                  typename internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type>
    {
    private:
        using Self = BasketDiff;
    public:
        using Base = typename internal::BasketDiffAggregate<BasketType,Type,Ext0,Exts...>::type;
        /// Weight function
        using WeightFunction = typename BasketType::WeightFunction;
        /// Point type used for computation
        using DataPoint = typename BasketType::DataPoint;
        /// Scalar type used for computation, as defined from Basket
        using Scalar = typename BasketType::Scalar;
    WRITE_COMPUTE_FUNCTIONS

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
    \copydoc BasketComputeObject
    \tparam P Implements \ref ponca_concepts "PointConcept"
    \tparam W Implements \ref concepts_weighting "WeightKernelConcept"
    \tparam Ext0 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept"
    \tparam Exts Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
*/
    template <class P, class W,
        template <class, class, typename> class Ext0,
        template <class, class, typename> class... Exts>
    class Basket : public BasketComputeObject<Basket<P, W, Ext0, Exts...>,
                                              typename internal::BasketAggregate<P, W, Ext0, Exts...>::type>
    {
    private:
        using Self = Basket;
    public:
        using Base = typename internal::BasketAggregate<P,W,Ext0,Exts...>::type;
        /// Weight function
        using WeightFunction = W;
        /// Point type used for computation
        using DataPoint = P;
        /// Scalar type used for computation, as defined from Basket
        using Scalar = typename DataPoint::Scalar;

        WRITE_COMPUTE_FUNCTIONS

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

#undef WRITE_COMPUTE_FUNCTIONS
} //namespace Ponca

