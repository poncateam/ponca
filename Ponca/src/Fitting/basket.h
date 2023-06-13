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

#ifndef PARSED_WITH_DOXYGEN
/*! \brief Namespace used for structure or classes used internally by the lib */
namespace internal
{
    /*! \brief Internal class used to build the Basket structure */
    template <class, class, typename T> class Forward: public T {};
    /*! \brief Internal class used to build the BasketDiff structure */
    template <class, class, int Type, typename T> class ForwardDiff: public T {};
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
    /*! The fit is evaluated multiple time if needed (see NEED_OTHER_PASS)*/                          \
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

#define BASKETDIFF_TP(I) template <class, class, int, typename> class Ext##I = internal::ForwardDiff
#define BSKW typename BasketType::WeightFunction
#define BSKP typename BasketType::DataPoint


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
         \tparam Ext1 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext2 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext3 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext4 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext5 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext6 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext7 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext8 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext9 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext10 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
         \tparam Ext11 Implements \ref concepts_computObjectBasketDiff "ComputationalDerivativesConcept" (optional)
     */
    template < typename BasketType, int Type,
            BASKETDIFF_TP(0), BASKETDIFF_TP(1), BASKETDIFF_TP(2), BASKETDIFF_TP(3), BASKETDIFF_TP(4), BASKETDIFF_TP(5), BASKETDIFF_TP(6), BASKETDIFF_TP(7), BASKETDIFF_TP(8), BASKETDIFF_TP(9), BASKETDIFF_TP(10), BASKETDIFF_TP(11) >
    class BasketDiff
        : public Ext11<BSKP,BSKW,Type, Ext10<BSKP,BSKW,Type, Ext9<BSKP,BSKW,Type, Ext8<BSKP,BSKW,Type, Ext7<BSKP,BSKW,Type, Ext6<BSKP,BSKW,Type, Ext5<BSKP,BSKW,Type, Ext4<BSKP,BSKW,Type, Ext3<BSKP,BSKW,Type, Ext2<BSKP,BSKW,Type, Ext1<BSKP,BSKW,Type,Ext0<BSKP,BSKW, Type, PrimitiveDer<BSKP,BSKW,Type, BasketType > > > > > > > > > > > > > {
    private:
        using Self   = BasketDiff;
    public:
    /// Base type, which aggregates all the computational objects using the CRTP
    using Base = Ext11 <BSKP, BSKW, Type,
            Ext10< BSKP, BSKW, Type,
            Ext9 < BSKP, BSKW, Type,
            Ext8 < BSKP, BSKW, Type,
            Ext7 < BSKP, BSKW, Type,
            Ext6 < BSKP, BSKW, Type,
            Ext5 < BSKP, BSKW, Type,
            Ext4 < BSKP, BSKW, Type,
            Ext3 < BSKP, BSKW, Type,
            Ext2 < BSKP, BSKW, Type,
            Ext1 < BSKP, BSKW, Type,
            Ext0 < BSKP, BSKW, Type,
            PrimitiveDer<BSKP,BSKW,Type,
            BasketType > > > > > > > > > > > > >;
        /// Weighting function
    using WeightFunction = BSKW;
    /// Point type used for computation
    using DataPoint = BSKP;
    /// Scalar type used for computation, as defined from Basket
    using Scalar = typename DataPoint::Scalar;

    WRITE_BASKET_FUNCTIONS

    /// Hide Basket::addNeighbor
    PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei) {
        // compute weight
        auto wres = Base::m_w.w(_nei.pos(), _nei);
        typename Base::ScalarArray dw;

        if (wres.first > Scalar(0.)) {
            Base::addLocalNeighbor(wres.first, wres.second, _nei, dw);
            Base::m_sumW += (wres.first);
            ++(Base::m_nbNeighbors);
            return true;
        }
        return false;
    }
};

#undef BASKETDIFF_TP


#define BASKET_TP(I) template <class, class, typename> class Ext##I = internal::Forward
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
    \tparam Ext1 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext2 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext3 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext4 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext5 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext6 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext7 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext8 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext9 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext10 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
    \tparam Ext11 Implements \ref concepts_computObjectBasket "ComputationalObjectConcept" (optional)
*/
    template < class P, class W,
        BASKET_TP(0), BASKET_TP(1), BASKET_TP(2), BASKET_TP(3), BASKET_TP(4), BASKET_TP(5), BASKET_TP(6), BASKET_TP(7), BASKET_TP(8), BASKET_TP(9), BASKET_TP(10), BASKET_TP(11) >
    class Basket
        : public Ext11<P,W, Ext10<P,W, Ext9<P,W, Ext8<P,W, Ext7<P,W, Ext6<P,W, Ext5<P,W, Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, PrimitiveBase<P,W> > > > > > > > > > > > >
    {
    private:
        using Self   = Basket;
    public:
        /// Base type, which aggregates all the computational objects using the CRTP
        using Base = Ext11<P,W, Ext10<P,W, Ext9<P,W, Ext8<P,W, Ext7<P,W, Ext6<P,W, Ext5<P,W, Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, PrimitiveBase<P,W> > > > > > > > > > > > >;
        /// Scalar type used for computation, as defined from template parameter `P`
        using Scalar = typename P::Scalar;
        /// Point type used for computation
        using DataPoint = P;
        /// Weighting function
        using WeightFunction = W;

        WRITE_BASKET_FUNCTIONS;

        /// \brief Add a neighbor to perform the fit
        /// When using by hand, don't forget to set m_nbNeighbors to 0 when performing a new pass (see PONCA::NEED_OTHER_PASS).
        /// \return false if param nei is not a valid neighbor (weight = 0)
        PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei) {
            // compute weight
            auto wres = Base::m_w.w(_nei.pos(), _nei);

            if (wres.first > Scalar(0.)) {
                Base::addLocalNeighbor(wres.first, wres.second, _nei);
                Base::m_sumW += (wres.first);
                ++(Base::m_nbNeighbors);
                return true;
            }
            return false;
        }
    }; // class Basket

#undef BASKET_TP

}// namespace Ponca

