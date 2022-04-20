/*
 Copyright (C) 2021 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include "./defines.h"
#include "./primitive.h"

#include <Eigen/Dense>

namespace Ponca
{

/*!
   \brief Empty fitting object doing no computation

   \inherit Concept::FittingProcedureConcept

   \ingroup fitting
 */

    template < class DataPoint, class _WFunctor, typename T>
    class DryFit :  public PrimitiveBase<DataPoint, _WFunctor>
    {
    private:
        typedef PrimitiveBase<DataPoint, _WFunctor> Base;

    public:
        /*! \brief Scalar type inherited from DataPoint*/
        typedef typename DataPoint::Scalar     Scalar;
        /*! \brief Vector type inherited from DataPoint*/
        typedef typename DataPoint::VectorType VectorType;
        /*! \brief Weight Function*/
        typedef _WFunctor                 WFunctor;

    public:
        /*! \brief Default constructor */
        PONCA_MULTIARCH inline DryFit() : Base() {}
        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
        PONCA_MULTIARCH inline void setWeightFunc (const WFunctor& /*_w*/) { }

        /*! \copydoc Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init (const VectorType& /*_evalPos*/) { }

        /**************************************************************************/
        /* Processing                                                             */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::addLocalNeighbor() */
        PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes)
        { return Base::addLocalNeighbor(w, localQ, attributes);}

        /*! \copydoc Concept::FittingProcedureConcept::finalize() */
        PONCA_MULTIARCH inline FIT_RESULT finalize() { return STABLE; }

        //! \brief Simulate Scalar field computation
        PONCA_MULTIARCH inline Scalar potential ( ) const { return Scalar(0); }

        //! \brief Simulate Scalar field computation
        PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const { return Scalar(0); }

        //! \brief Simulate point projection
        PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const { return _q; }

        //! \brief Simulate gradient direction computation
        PONCA_MULTIARCH inline VectorType primitiveGradient () const { return VectorType::Zero(); }

        //! \brief Simulate gradient direction computation
        PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType&) const { return VectorType::Zero(); }
    };

} //namespace Ponca
