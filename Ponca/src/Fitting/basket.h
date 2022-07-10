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

/*! \brief Namespace used for structure or classes used internally by the lib */
namespace internal
{
    /*! \brief Internal class used to build the Basket structure */
    template <class, class, typename T> class Forward: public T {};
}


#define BASKET_TP(I) template <class, class, typename> class Ext##I = internal::Forward
/*!

    \brief Agregator class used to declare specialized structures using CRTP (Curiously Recurring Template Pattern)
    \todo Comment

    The various implementations of Ponca::Concept are mixed through
    specializations of the Basket class:
    \code
        typedef
        Basket <PointImpl,      // Implementation of PointConcept
        WeightFuncImpl,         // Implementation of WeightFuncConcept (could use WeightKernelConcept)
        FittingProcedureImpl,   // Implementation of FittingProcedureConcept
        FittingExtensionImpl1,  //
        FittingExtensionImpl2,  // Implementations of FittingExtensionConcept
        ... ,                   //
        FittingExtensionImpln   //
        > myFit;                // Final structure to fit and extend a primitive over weighted samples
    \endcode

    \ingroup fitting
*/
    template < class P, class W, template <class, class, typename>class Fit,
        BASKET_TP(0), BASKET_TP(1), BASKET_TP(2), BASKET_TP(3), BASKET_TP(4), BASKET_TP(5), BASKET_TP(6), BASKET_TP(7), BASKET_TP(8), BASKET_TP(9), BASKET_TP(10), BASKET_TP(11) >
    class Basket
        : public Ext11<P,W, Ext10<P,W, Ext9<P,W, Ext8<P,W, Ext7<P,W, Ext6<P,W, Ext5<P,W, Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, Fit<P,W, PrimitiveBase<P,W> > > > > > > > > > > > > >
    {
    public:
        using Base = Ext11<P,W, Ext10<P,W, Ext9<P,W, Ext8<P,W, Ext7<P,W, Ext6<P,W, Ext5<P,W, Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, Fit<P,W, PrimitiveBase<P,W> > > > > > > > > > > > > >;
        using Scalar = typename P::Scalar;
        typedef P DataPoint;
        typedef W WeightFunction;
        /*!
         * \brief Convenience function for STL-like iterators
         *
         * Add neighbors stored in a container using STL-like iterators, and
         * call finalize at the end.
         * The fit is evaluated multiple time if needed (see NEED_OTHER_PASS).
         */
        template <typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline
        FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end){
            FIT_RESULT res = UNDEFINED;
            do {
                for (auto it = begin; it != end; ++it){
                    addNeighbor(*it);
                }
                res = Base::finalize();
            } while ( res == NEED_OTHER_PASS );
            return res;
        }



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

        /*!
         * \brief Convenience function for STL-like iterators
         *
         * Add neighbors stored in a PointContainer and sampled using
         * indices stored in ids.
         *
         * \tparam IndexRange STL-Like range storing indices of the neighbors
         * \tparam PointContainer STL-like container storing the points
         *
         * \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
         */
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline
        FIT_RESULT computeWithIds(IndexRange ids, const PointContainer& points){
            FIT_RESULT res = UNDEFINED;
            do {
                for (const auto& i : ids){
                    this->addNeighbor(points[i]);
                }
                res = this->finalize();
            } while ( res == NEED_OTHER_PASS );
            return res;
        }

#ifdef PONCA_CPU_ARCH
        /*!
         * \brief Convenience function for STL-like containers
         *
         * \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
         */
        template <typename Container>
        PONCA_MULTIARCH inline
        FIT_RESULT compute(const Container& c){
            return compute(std::begin(c), std::end(c));
        }
#endif
    }; // class Basket

#undef BASKET_TP

}// namespace Ponca

