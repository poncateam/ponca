/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>


 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./surfacePrimitive.h"
#include <Eigen/Eigenvalues>


#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

#include <Eigen/Core>

namespace Ponca
{

    /*!
        \brief Surface fitting procedure using only points position

        This class can also computes the surface variation measure introduced in
        \cite Pauly:2002:PSSimplification. The solver used to analyse the covariance
        matrix is stored for further use.

        \inherit Concept::FittingProcedureConcept

        \see CompactSurface
        \ingroup fitting
    */
    template < class DataPoint, class _WFunctor, typename T >
    class LeastSquareSurfaceFit : public Surface<DataPoint, _WFunctor>
    {
    private:
        typedef Surface<DataPoint, _WFunctor> Base;

    protected:
        enum
        {
            /*!
            * \brief Expose a method worldToTangentSurface(VectorType), which turns a point
            * in ambient 3D space to the tangent surface.
            * \see worldToTangentSurface
            * \see tangentSurfaceToWorld
            */
            PROVIDES_TANGENT_SURFACE_BASIS
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
        VectorType m_cog;    /*!< \brief Gravity center of the neighborhood */
        Eigen::Matrix<Scalar, 9, 1 > cofficient;   /*!< \brief cofficients of the fitting equation */


        Eigen::Matrix<Scalar, 9, 1>  m_right;   
        Eigen::Matrix<Scalar, 9, 9>  m_cov;     /*!< \brief Covariance matrix */

        Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */
        WFunctor m_w;     /*!< \brief Weight function (must inherits BaseWeightFunc) */

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline LeastSquareSurfaceFit() : Base() {}

        /*! \brief Explicit conversion to LeastSquareSurfaceFit, to access methods potentially hidden by inheritage */
        PONCA_MULTIARCH inline
        LeastSquareSurfaceFit<DataPoint, WFunctor, T>& leastSquareSurfaceFit()
        { return * static_cast<LeastSquareSurfaceFit<DataPoint, WFunctor, T>*>(this); }    

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

        /**************************************************************************/
        /* Results                                                                */
        /**************************************************************************/

        /*! \brief Reading access to the Solver used to analyse the covariance
        matrix */
        PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }

        /*!
        \brief Project a point on the surface using Gradient Descent
        This projection is realized by following the gradient of the surface scalar field
        \param nbIter Number of iterations (default = 16)
        */
        PONCA_MULTIARCH inline VectorType project (const VectorType& _q, int nbIter = 16) const;


        //! \brief Approximation of the scalar field gradient at \f$ \mathbf{q} (not normalized) \f$
        PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType& _q) const;

    }; 

    #include "leastSquareSurfaceFit.hpp"

} //namespace Ponca
