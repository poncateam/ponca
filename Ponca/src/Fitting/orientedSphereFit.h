/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./algebraicSphere.h"
#include "./mean.h"          // used to define OrientedSphereFit

namespace Ponca
{

/*!
    \brief Algebraic Sphere fitting procedure on oriented point sets

    Method published in \cite Guennebaud:2007:APSS.

    \inherit Concept::FittingProcedureConcept

    \see AlgebraicSphere
*/
template < class DataPoint, class _WFunctor, typename T >
class OrientedSphereFitImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE &&
                Base::PROVIDES_MEAN_NORMAL &&
                Base::PROVIDES_MEAN_POSITION
    };

    // computation data
    Scalar  m_sumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals */
            m_sumDotPP, /*!< \brief Sum of the squared relative positions */
            m_nume,     /*!< \brief Numerator of the quadratic parameter (excluding the 0.5 coefficient)   */
            m_deno;     /*!< \brief Denominator of the quadratic parameter (excluding the 0.5 coefficient) */

public:
    PONCA_EXPLICIT_CAST_OPERATORS(OrientedSphereFitImpl,orientedSphereFit)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
    PONCA_FITTING_IS_SIGNED(true)
}; //class OrientedSphereFitImpl

/// \brief Helper alias for Oriented Sphere fitting on 3D points using OrientedSphereFitImpl
//! [OrientedSphereFit Definition]
template < class DataPoint, class _WFunctor, typename T>
using OrientedSphereFit =
OrientedSphereFitImpl<DataPoint, _WFunctor,
        MeanPosition<DataPoint, _WFunctor,
            MeanNormal<DataPoint, _WFunctor,
                AlgebraicSphere<DataPoint, _WFunctor,T>>>>;
//! [OrientedSphereFit Definition]

/*!
    \brief Internal generic class performing the Fit derivation
*/
template < class DataPoint, class _WFunctor, int DiffType, typename T>
class OrientedSphereDerImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE &
                Base::PROVIDES_MEAN_POSITION_DERIVATIVE &
                Base::PROVIDES_PRIMITIVE_DERIVATIVE,
        PROVIDES_ALGEBRAIC_SPHERE_DERIVATIVE,
        PROVIDES_NORMAL_DERIVATIVE
    };

    // computation data
    VectorArray m_dSumN;     /*!< \brief Sum of the normal vectors with differenciated weights */
    ScalarArray m_dSumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals with differenciated weights */
                m_dSumDotPP, /*!< \brief Sum of the squared relative positions with differenciated weights */
                m_dNume,     /*!< \brief Differenciation of the numerator of the quadratic parameter   */
                m_dDeno;     /*!< \brief Differenciation of the denominator of the quadratic parameter */

public:
    // results
    ScalarArray m_dUc, /*!< \brief Derivative of the hyper-sphere constant term  */
                m_dUq; /*!< \brief Derivative of the hyper-sphere quadratic term */
    VectorArray m_dUl; /*!< \brief Derivative of the hyper-sphere linear term    */

public:
    PONCA_EXPLICIT_CAST_OPERATORS_DER(OrientedSphereDerImpl,orientedSphereDer)
    PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE

    /*! \brief Returns the derivatives of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline ScalarArray dPotential() const;

    /*! \brief Returns the derivatives of the primitive normal */
    PONCA_MULTIARCH inline VectorArray dNormal() const;

    /*! \brief compute  the square of the Pratt norm derivative */
    PONCA_MULTIARCH inline ScalarArray dprattNorm2() const
    {
        return   Scalar(2.) * Base::m_ul.transpose() * m_dUl
            - Scalar(4.) * Base::m_uq * m_dUc
            - Scalar(4.) * Base::m_uc * m_dUq;
    }

    /*! \brief compute the square of the Pratt norm derivative for dimension _d */
    PONCA_MULTIARCH inline Scalar dprattNorm2(unsigned int _d) const
    {
        return   Scalar(2.) * m_dUl.col(_d).dot(Base::m_ul)
            - Scalar(4.) * m_dUc.col(_d)[0]*Base::m_uq
            - Scalar(4.) * m_dUq.col(_d)[0]*Base::m_uc;
    }

    /*! \brief compute the Pratt norm derivative for the dimension _d */
    PONCA_MULTIARCH inline Scalar dprattNorm(unsigned int _d) const
    {
        PONCA_MULTIARCH_STD_MATH(sqrt);
        return sqrt(dprattNorm2(_d));
    }

    /*! \brief compute the Pratt norm derivative */
    PONCA_MULTIARCH inline Scalar dprattNorm() const
    {
        PONCA_MULTIARCH_STD_MATH(sqrt);
        return dprattNorm2().array().sqrt();
    }
    //! Normalize the scalar field by the Pratt norm
    /*!
        \warning Requieres that isNormalized() return false
        \return false when the original sphere has already been normalized.
    */
    PONCA_MULTIARCH inline bool applyPrattNorm();

}; //class OrientedSphereDerImpl

/// \brief Helper alias for Oriented Sphere fitting on 3D points using OrientedSphereDerImpl
    template < class DataPoint, class _WFunctor, int DiffType, typename T>
    using OrientedSphereDer =
        OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType,
            MeanPositionDer<DataPoint, _WFunctor, DiffType, T>>;

#include "orientedSphereFit.hpp"

} //namespace Ponca
