/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./algebraicSphere.h"


namespace Ponca
{

/*!
    \brief Algebraic Sphere fitting procedure on point sets without normals

    Method published in \cite Guennebaud:2007:APSS.

    \inherit Concept::FittingProcedureConcept

    \see AlgebraicSphere

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T >
class SphereFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE
    };
protected:
    typedef Eigen::Matrix<Scalar, DataPoint::Dim+2, 1>      VectorA;
    typedef Eigen::Matrix<Scalar, DataPoint::Dim+2, DataPoint::Dim+2>  MatrixA;

    // computation data
    MatrixA  m_matA {MatrixA::Zero()};  /*!< \brief Covariance matrix of [1, p, p^2] */

public:
    PONCA_EXPLICIT_CAST_OPERATORS(SphereFitImpl,sphereFit)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
}; //class SphereFit

/// \brief Helper alias for Sphere fitting on 3D points using SphereFitImpl
/// \ingroup fittingalias
template < class DataPoint, class _WFunctor, typename T>
using SphereFit =
    SphereFitImpl<DataPoint, _WFunctor,
        AlgebraicSphere<DataPoint, _WFunctor,T>>;



#include "sphereFit.hpp"

} //namespace Ponca
