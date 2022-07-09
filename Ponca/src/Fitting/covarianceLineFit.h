/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>

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
   \brief Line fitting procedure that minimize the orthogonal distance between the samples and the fitted primitive.

   \inherit Concept::FittingProcedureConcept
   \see Line
   \see CovariancePlaneFit which use a similar approach for Plane estimation

   \warning This class is valid only in 3D.
   \ingroup fitting
 */

template < class DataPoint, class _WFunctor, typename T>
class CovarianceLineFit : public Line<DataPoint, _WFunctor>
{
private:
    typedef Line<DataPoint, _WFunctor> Base;

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using MatrixType = typename Base::MatrixType; /*!< \brief Inherited matrix type*/
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/
    /*! \brief Solver used to analyse the covariance matrix*/
    typedef Eigen::SelfAdjointEigenSolver<MatrixType> Solver;

protected:
     // computation data
    VectorType m_cog {VectorType::Zero()};     /*!< \brief Gravity center of the neighborhood */
    MatrixType m_cov {MatrixType::Zero()};     /*!< \brief Covariance matrix */

    Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */

public:
     /*! \brief Default constructor */
    PONCA_MULTIARCH inline CovarianceLineFit() = default;
    //! \brief Explicit conversion to CovarianceLineFit, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    CovarianceLineFit<DataPoint, WFunctor, T>& leastSquareLine()
    { return * static_cast<CovarianceLineFit<DataPoint, WFunctor, T>*>(this); }
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

    /*! \brief Reading access to the Solver used to analyse the covariance
      matrix */
    PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }
  

};

#include "covarianceLineFit.hpp"

} //namespace Ponca
