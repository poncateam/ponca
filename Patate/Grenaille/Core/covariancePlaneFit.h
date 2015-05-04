/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_COVARIANCE_PLANE_FIT_
#define _GRENAILLE_COVARIANCE_PLANE_FIT_

#include <Eigen/Eigenvalues> 

namespace Grenaille
{

/*!
    \brief Plane fitting procedure using only points position
    
    \note This procedure requires two passes to fit a plane

    This class can also computes the surface variation measure introduced in
    \cite Pauly:2002:PSSimplification. The solver used to analyse the covariance 
    matrix is stored for further use.

    \inherit Concept::FittingProcedureConcept

    \warning This class is currently untested and should not be used !

    \see CompactPlane
*/
template < class DataPoint, class _WFunctor, typename T >
class CovariancePlaneFit : public T
{
private:
    typedef T Base;

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE
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
    Scalar  m_sumW;       /*!< \brief Sum of queries weight.*/
    VectorType m_gc,      /*!< \brief Gravity center of the neighborhood */
               m_evalPos; /*!< \brief Center of the evaluation basis */
    MatrixType m_cov;     /*!< \brief Covariance matrix */
    /*! \brief State indicating if the current pass is the first one*/
    bool       m_isFirstPass; 

    Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */

    WFunctor m_w;     /*!< \brief Weight function (must inherits BaseWeightFunc) */

public:

    /*! \brief Default constructor */
    MULTIARCH inline CovariancePlaneFit() : Base() {}

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
    MULTIARCH inline void setWeightFunc (const WFunctor& _w) { m_w  = _w; }

    /*! \copydoc Concept::FittingProcedureConcept::init() */
    MULTIARCH inline void init (const VectorType& _evalPos);

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
    MULTIARCH inline bool addNeighbor(const DataPoint &_nei);

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    MULTIARCH inline FIT_RESULT finalize();

    /**************************************************************************/
    /* Results                                                                */
    /**************************************************************************/
    /*! \brief Reading access to the Solver used to analyse the covariance 
      matrix */
    MULTIARCH inline const Solver& solver() const { return m_solver; }
    
    /*! \brief Implements \cite Pauly:2002:PSSimplification. 
        \return 0 for invalid fits */
    MULTIARCH inline Scalar surfaceVariation() const;
}; //class OrientedSphereFit


#include "covariancePlaneFit.hpp"

} //namespace Grenaille


#endif
