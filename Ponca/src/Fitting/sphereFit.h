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
private:
    using Base = T;

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/

protected:
    typedef Eigen::Matrix<Scalar, DataPoint::Dim+2, 1>      VectorA;
    typedef Eigen::Matrix<Scalar, DataPoint::Dim+2, DataPoint::Dim+2>  MatrixA;

    // computation data
    MatrixA  m_matA {MatrixA::Zero()};  /*!< \brief Covariance matrix of [1, p, p^2] */

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline SphereFitImpl() = default;

    PONCA_EXPLICIT_CAST_OPERATORS(SphereFitImpl,sphereFit)

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/

    /*! \copydoc Concept::FittingProcedureConcept::init() */
    PONCA_MULTIARCH inline void init (const VectorType& _evalPos);


    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addLocalNeighbor() */
    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize();

}; //class SphereFit

/// \brief Helper alias for Sphere fitting on 3D points using SphereFitImpl
/// \ingroup fittingalias
template < class DataPoint, class _WFunctor, typename T>
using SphereFit =
    SphereFitImpl<DataPoint, _WFunctor,
        AlgebraicSphere<DataPoint, _WFunctor,T>>;



#include "sphereFit.hpp"

} //namespace Ponca
