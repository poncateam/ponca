/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 Copyright (C) 2022 nicolas mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./mean.h"
#include <Eigen/Dense>

namespace Ponca
{

/*!
   \brief Procedure that compute and decompose the covariance matrix of the neighbors positions in \f$3d\f$.

   This process is commonly used for plane fitting and local variance analysis. It is often called Principal
   Component Analysis (PCA) of the neighborhood, and used in Geometry Processing and Computer Vision.

   \inherit Concept::FittingProcedureConcept
   \see CovariancePlaneFit which use a similar approach for Plane estimation

   ### Computation details
   Standard PCA algorithm involves a two-steps process where the barycenter \f$\mathbf{b}\f$ is first computed,
   and then the covariance matrix \f$\text{C}\f$ (in the following, the weights are ignored for clarity but without
   loss of generality):
   \f{align}
   \mathbf{b} &= \frac{1}{n}\sum_i\mathbf{p}_i \\
   \text{C}   &= \frac{1}{n}\sum_i(\mathbf{p}_i-\mathbf{b})(\mathbf{p}_i-\mathbf{b})^T
   \f}
   This class implements a single-pass version, where the first formulation is re-expressed as follows:
   \f{align}
   \text{C} &= \frac{1}{n}\sum_i (\mathbf{p}_i\mathbf{p}_i^T - \mathbf{b}\mathbf{p}_i^T - \mathbf{p}_i\mathbf{b}^T + \mathbf{b}\mathbf{b}^T) \\
            &= \frac{1}{n}\sum_i (\mathbf{p}_i\mathbf{p}_i^T) -  \frac{1}{n}\sum_i(\mathbf{b}\mathbf{p}_i^T)  -  \frac{1}{n}\sum_i(\mathbf{p}_i\mathbf{b}^T)  +  \frac{1}{n}\sum_i (\mathbf{b}\mathbf{b}^T) \\
            &= \frac{1}{n}\sum_i (\mathbf{p}_i\mathbf{p}_i^T) - \mathbf{b}\frac{1}{n}\sum_i(\mathbf{p}_i^T) - \frac{1}{n}\sum_i(\mathbf{p}_i)\mathbf{b}^T  + \frac{1}{n}\sum_i(1) \mathbf{b}\mathbf{b}^T \\
            &= \frac{1}{n}\sum_i (\mathbf{p}_i\mathbf{p}_i^T) - \mathbf{b}\mathbf{b}^T - \mathbf{b}\mathbf{b}^T + \mathbf{b}\mathbf{b}^T \f}
   Leading to a single pass where \f$\text{C}\f$ is express by substracting two terms that can be computed independently
   in one run:
   \f[ \text{C} = \frac{1}{n}\sum_i (\mathbf{p}_i\mathbf{p}_i^T) - \mathbf{b}\mathbf{b}^T \f]

   All the computed features are defined for the 3 eigenvalues \f$ 0 < \lambda_0
   \leq \lambda_1 \leq \lambda_2 \f$.


   \warning This class is valid only in 3D.
 */

    template < class DataPoint, class _WFunctor, typename T>
    requires ProvidesMeanPosition<T>
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

        /*! \brief Implements the planarity \cite Guinard:2017 .
            Planarity is defined as:
            \f[ \frac{\lambda_1 - \lambda_0}{\lambda_2} \f]
        */
        PONCA_MULTIARCH inline Scalar planarity() const;

        /*! \brief Implements the linearity \cite Guinard:2017 .
            Linearity is defined as:
            \f[ \frac{\lambda_2 - \lambda_1}{\lambda_2} \f]
        */
        PONCA_MULTIARCH inline Scalar linearity() const;

        /*! \brief Implements the sphericity \cite Guinard:2017 .
            Sphericity is defined as:
            \f[ \frac{\lambda_0}{\lambda_2} \f]
        */
        PONCA_MULTIARCH inline Scalar sphericity() const;

        /*! \brief Implements the anisotropy \cite Guinard:2017 .
            Anisotropy is defined as:
            \f[ \frac{\lambda_2 - \lambda_0}{\lambda_2} \f]
        */
        PONCA_MULTIARCH inline Scalar anisotropy() const;

        /*! \brief Implements the eigenentropy \cite Guinard:2017 .
            Eigenentropy is defined as:
            \f[ - \lambda_0 * \ln{\lambda_0} - \lambda_1 * \ln{\lambda_1} - \lambda_2 * \ln{\lambda_2} \f]
        */
        PONCA_MULTIARCH inline Scalar eigenentropy() const;

        /*! \brief The minimun eigenvalue \f$ \lambda_0 \f$.
        */
        PONCA_MULTIARCH inline Scalar lambda_0() const;

        /*! \brief The second eigenvalue \f$ \lambda_1 \f$.
        */
        PONCA_MULTIARCH inline Scalar lambda_1() const;

        /*! \brief The maximun eigenvalue \f$ \lambda_2 \f$.
        */
        PONCA_MULTIARCH inline Scalar lambda_2() const;

        /*! \brief Reading access to the Solver used to analyse the covariance matrix */
        PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }
    };


/*!
    \brief Internal generic class computing the derivatives of covariance matrix
    computed by CovarianceFitBase
    \inherit Concept::FittingExtensionConcept
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
