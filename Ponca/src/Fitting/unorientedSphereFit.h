/*
 Copyright (C) 2013 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./algebraicSphere.h"
#include "./mean.h"          // used to define UnorientedSphereFit

#include <Eigen/Dense>

namespace Ponca
{
/*!
    \brief Algebraic Sphere fitting procedure on point sets with non-oriented normals

    Method published in \cite Chen:2013:NOMG.

    \inherit Concept::FittingProcedureConcept

    \see class AlgebraicSphere, class OrientedSphereFit
*/
template < class DataPoint, class _WFunctor, typename T >
class UnorientedSphereFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum { Check = Base::PROVIDES_ALGEBRAIC_SPHERE && Base::PROVIDES_MEAN_POSITION };

    typedef Eigen::Matrix<Scalar, DataPoint::Dim+1, 1>      VectorB;
    typedef Eigen::Matrix<Scalar, DataPoint::Dim+1, DataPoint::Dim+1>  MatrixBB;
    
public:
    using Solver = Eigen::EigenSolver<MatrixBB>;

    MatrixBB    m_matA {MatrixBB::Zero()}; /*!< \brief The accumulated covariance matrix */
    MatrixBB    m_matQ {MatrixBB::Zero()}; /*!< \brief The constraint matrix */
    Scalar      m_sumDotPP {0};            /*!< \brief Sum of the squared relative positions */

    Solver m_solver;
    
public:
    PONCA_EXPLICIT_CAST_OPERATORS(UnorientedSphereFitImpl,unorientedSphereFit)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
    PONCA_FITTING_IS_SIGNED(false)

}; // class UnorientedSphereFitImpl

/// \brief Helper alias for Oriented Sphere fitting on 3D points using UnorientedSphereFitImpl
template < class DataPoint, class _WFunctor, typename T>
using UnorientedSphereFit =
UnorientedSphereFitImpl<DataPoint, _WFunctor,
        MeanPosition<DataPoint, _WFunctor,
                AlgebraicSphere<DataPoint, _WFunctor,T>>>;

} //namespace Ponca

#include "unorientedSphereFit.hpp"
