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

    This method published in \cite Chen:2013:NOMG maximizes the sum of squared dot product between the input normal vectors and the gradient of the algebraic sphere.
    The maximization is done under the constraint that the norm of the gradient is unitary on average.
    In practice, it amounts to solve the generalized eigenvalue problem
    \f[
        A \mathbf{u} = \lambda Q \mathbf{u}
    \f]
    where
    \f[
        \mathbf{u} = \begin{bmatrix} u_l \\ u_q \end{bmatrix}
    \f]
    \f[
        A = \sum_i w_i \begin{bmatrix} \mathbf{n}_i \\ \mathbf{n}_i^T\mathbf{p}_i \end{bmatrix}\begin{bmatrix} \mathbf{n}_i^T & \mathbf{n}_i^T\mathbf{p}_i \end{bmatrix}
    \f]
    \f[
        Q = \frac{1}{\sum_i w_i} \begin{bmatrix} I & \sum_i w_i \mathbf{p}_i \\ \sum_i w_i \mathbf{p}_i^T & \sum_i w_i \mathbf{p}_i^T\mathbf{p}_i \end{bmatrix}
    \f]
    The constant coefficient \f$u_c\f$ is computed as in \ref OrientedSphereFitImpl by minimizing the sum of squared potential evaluated at the points \f$\mathbf{p}_i\f$.
    This fitting method corresponds to the case where the scalar field is \f$ f^* \f$ in the Section 3.3 of \cite Chen:2013:NOMG.

    \inherit Concept::FittingProcedureConcept

    \see class AlgebraicSphere, class OrientedSphereFit
*/
template < class DataPoint, class _NFilter, typename T >
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
template < class DataPoint, class _NFilter, typename T>
using UnorientedSphereFit =
UnorientedSphereFitImpl<DataPoint, _NFilter,
        MeanPosition<DataPoint, _NFilter,
                AlgebraicSphere<DataPoint, _NFilter, T>>>;



template < class DataPoint, class _NFilter, int DiffType, typename T>
class UnorientedSphereDerImpl : public T
{
protected:
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

    using VectorB = typename Base::VectorB;
    using MatrixBB = typename Base::MatrixBB;

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE &
                Base::PROVIDES_MEAN_POSITION_DERIVATIVE &
                Base::PROVIDES_PRIMITIVE_DERIVATIVE,
        PROVIDES_ALGEBRAIC_SPHERE_DERIVATIVE,
        PROVIDES_NORMAL_DERIVATIVE
    };

protected:
    // computation data
    MatrixBB m_dmatA[Base::NbDerivatives];
    ScalarArray m_dSumDotPP;

public:
    // results
    ScalarArray m_dUc;
    VectorArray m_dUl;
    ScalarArray m_dUq;

public:
    PONCA_EXPLICIT_CAST_OPERATORS_DER(UnorientedSphereDerImpl,unorientedSphereDer)
    PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE

    PONCA_MULTIARCH inline ScalarArray dPotential() const;
    PONCA_MULTIARCH inline VectorArray dNormal() const;

}; //class UnorientedSphereDerImpl


template < class DataPoint, class _NFilter, int DiffType, typename T>
using UnorientedSphereDer =
    UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType,
        MeanPositionDer<DataPoint, _NFilter, DiffType, T>>;


} //namespace Ponca

#include "unorientedSphereFit.hpp"
