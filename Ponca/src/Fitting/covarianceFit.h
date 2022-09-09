/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 Copyright (C) 2022 nicolas mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include "./defines.h"

#include <Eigen/Dense>

namespace Ponca
{

/*!
   \brief Line fitting procedure that minimize the orthogonal distance
   between the samples and the fitted primitive.

   \inherit Concept::FittingProcedureConcept
   \see Line
   \see CovariancePlaneFit which use a similar approach for Plane estimation

   \todo Add equations

   \warning This class is valid only in 3D.
   \ingroup fitting
 */

    template < class DataPoint, class _WFunctor, typename T>
    class CovarianceFitBase : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

    protected:
        enum
        {
            Check = Base::PROVIDES_MEAN_POSITION,
            PROVIDES_POSITION_COVARIANCE
        };

    public:
        using MatrixType = typename DataPoint::MatrixType; /*!< \brief Alias to matrix type*/
        /*! \brief Solver used to analyse the covariance matrix*/
        using Solver = Eigen::SelfAdjointEigenSolver<MatrixType>;

    protected:
        // computation data
        MatrixType m_cov {MatrixType::Zero()};     /*!< \brief Covariance matrix */
        Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(CovarianceFitBase,covarianceFit)
        PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

        /*! \brief Implements \cite Pauly:2002:PSSimplification surface variation.
            It computes the ratio \f$ d \frac{\lambda_0}{\sum_i \lambda_i} \f$ with \c d the dimension of the ambient space.
            \return 0 for invalid fits
        */
        PONCA_MULTIARCH inline Scalar surfaceVariation() const;

        /*! \brief Reading access to the Solver used to analyse the covariance matrix */
        PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }
    };


/*!
    \brief Internal generic class computing the derivatives of covariance matrix
    computed by CovarianceFitBase
    \inherit Concept::FittingExtensionConcept

    \ingroup fitting
*/
    template < class DataPoint, class _WFunctor, int DiffType, typename T>
    class CovarianceFitDer : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE
    PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

    protected:
        enum
        {
            Check = Base::PROVIDES_PRIMITIVE_DERIVATIVE &&
                    Base::PROVIDES_MEAN_POSITION_DERIVATIVE &&
                    Base::PROVIDES_POSITION_COVARIANCE,
            PROVIDES_POSITION_COVARIANCE_DERIVATIVE
        };

    protected:
        /// Computation data: derivatives of the covariance matrix
        MatrixType  m_dCov[Base::NbDerivatives];

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(CovarianceFitDer,covarianceFitDer)
        PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE
    }; //class CovarianceFitDer

#include "covarianceFit.hpp"

} //namespace Ponca
