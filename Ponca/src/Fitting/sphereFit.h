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
    \brief Algebraic Sphere fitting procedure on point set without normals

    This method published in \cite Guennebaud:2007:APSS minimizes
    \f[
        \mathcal{L}(\mathbf{u}) = \frac{1}{2} \sum_i w_i f_{\mathbf{u}}(\mathbf{x}_i)^2 = \frac{1}{2} \mathbf{u}^T A \mathbf{u}
    \f]
    with \f$ A = \sum_i w_i \tilde{\mathbf{x}_i}  \tilde{\mathbf{x}_i}^T\f$ and \f$ f_{\mathbf{u}} \f$ the algebraic sphere defined by
    \f[
        f_{\mathbf{u}}(\mathbf{x}) =
        u_c + \mathbf{u}_l.\mathbf{x} + u_q \mathbf{x}.\mathbf{x} =
        \begin{bmatrix}
        1 & \mathbf{x}^T & \mathbf{x}.\mathbf{x}
        \end{bmatrix}
        \begin{bmatrix}
        u_c \\ u_l \\ u_q
        \end{bmatrix}
        = \tilde{\mathbf{x}}^T \mathbf{u},
    \f]
    under the constraint (unitary gradient onto the surface)
    \f[
        \|\mathbf{u}_l\|^2 - 4 u_c u_q = \mathbf{u}^T C \mathbf{u} = 1
    \f]
    where
    \f[
        C =
        \begin{bmatrix}
            0 &   &         &   & -2 \\
              & 1 &         &   &    \\
              &   &  \ddots &   &    \\
              &   &         & 1 &    \\
           -2 &   &         &   &  0
        \end{bmatrix}
    \f]
    which amounts to solve the generalized eigenvalue problem \f$ A\mathbf{u} = \lambda C \mathbf{u} \f$.

    \inherit Concept::FittingProcedureConcept

    \see AlgebraicSphere
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

public:
    using Solver = Eigen::EigenSolver<MatrixA>;

protected:
    // computation data
    MatrixA  m_matA {MatrixA::Zero()};  /*!< \brief Covariance matrix of [1, p, p^2] */

    Solver m_solver;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(SphereFitImpl,sphereFit)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }

}; //class SphereFit

/// \brief Helper alias for Sphere fitting on 3D points using SphereFitImpl
template < class DataPoint, class _WFunctor, typename T>
using SphereFit =
    SphereFitImpl<DataPoint, _WFunctor,
        AlgebraicSphere<DataPoint, _WFunctor,T>>;



#include "sphereFit.hpp"

} //namespace Ponca
