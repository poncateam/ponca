/*
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
 * \brief Extension to compute line primitives and fitting procedures to the library
 *
 * \ingroup fitting
 */

template < class DataPoint, class _WFunctor, typename T>
class LeastSquareLine : public Line<DataPoint, _WFunctor>
{
private:
    typedef Line<DataPoint, _WFunctor> Base;


protected:
    enum
    {
        Check = Base::PROVIDES_LINE
    };

public:
   /*! \brief Scalar type inherited from DataPoint*/
    typedef typename Base::Scalar     Scalar;
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename Base::VectorType VectorType;
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename Base::MatrixType MatrixType;
    /*! \brief Weight Function*/
    typedef _WFunctor                 WFunctor;
    /*! \brief Solver used to analyse the covariance matrix*/
    typedef Eigen::SelfAdjointEigenSolver<MatrixType> Solver;

protected:
     // computation data
    Scalar  m_sum;       /*!< \brief total number of queries .*/
    VectorType m_cog;     /*!< \brief Gravity center of the neighborhood */
    MatrixType m_cov;     /*!< \brief Covariance matrix */

    Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */
    WFunctor m_w;     /*!< \brief Weight function (must inherits BaseWeightFunc) */

public:
     /*! \brief Default constructor */
    PONCA_MULTIARCH inline LeastSquareLine() : Base() {}
    //! \brief Explicit conversion to LeastSquareLine, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    LeastSquareLine<DataPoint, WFunctor, T>& leastSquareLine()
    { return * static_cast<LeastSquareLine<DataPoint, WFunctor, T>*>(this); }
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
    PONCA_MULTIARCH inline void setWeightFunc (const WFunctor& _w) { m_w  = _w; }

    /*! \copydoc Concept::FittingProcedureConcept::init() */
    PONCA_MULTIARCH inline void init (const VectorType& _evalPos);

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
    PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei);

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize();

    /*! \brief Reading access to the Solver used to analyse the covariance
      matrix */
    PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }
  

};

#include "LeastSquareLine.hpp"

} //namespace Ponca
