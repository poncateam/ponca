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
        Check = Base::PROVIDES_SURFACE,
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

    /*! \brief Implements \cite Pauly:2002:PSSimplification surface variation.

        It computes the ratio \f$ d \frac{\lambda_0}{\sum_i \lambda_i} \f$ with \c d the dimension of the ambient space.

        \return 0 for invalid fits
    */
    PONCA_MULTIARCH inline Scalar surfaceVariation() const;

    /*!
     * \brief Express a point in ambiant space relatively to the tangent surface.
     * Output vector is: [h, u, v]^T, where u, v are 2d coordinates on the surface,
     * and h the height of the sample.
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType worldToTangentSurface(const VectorType &_q) const;

    /*!
     * \brief Transform a point from the tangent surface [h, u, v]^T to ambiant space
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType tangentSurfaceToWorld(const VectorType &_q) const;
}; //class LeastSquareSurfaceFit

// namespace internal {

// using ::Ponca::internal::FitSpaceDer;
// using ::Ponca::internal::FitScaleDer;

// /*!
//     \brief Internal generic class computing the derivatives of covariance surface fits
//     \inherit Concept::FittingExtensionConcept

//     The differentiation can be done automatically in scale and/or space, by
//     combining the enum values FitScaleDer and FitSpaceDer in the template
//     parameter Type.

//     The differenciated values are stored in static arrays. The size of the
//     arrays is computed with respect to the derivation type (scale and/or space)
//     and the number of the dimension of the ambiant space.
//     By convention, the scale derivatives are stored at index 0 when Type
//     contains at least FitScaleDer. The size of these arrays can be known using
//     derDimension(), and the differentiation type by isScaleDer() and
//     isSpaceDer().

//     \ingroup fitting
//     \internal
// */
// template < class DataPoint, class _WFunctor, typename T, int Type>
// class CovarianceSurfaceDer : public T
// {
// private:
//     typedef T Base; /*!< \brief Generic base type */


// protected:
//     enum
//     {
//         Check = Base::PROVIDES_SURFACE,             /*!< \brief Needs surface */
//         PROVIDES_COVARIANCE_SURFACE_DERIVATIVE,      /*!< \brief Provides derivatives for hyper-surfaces */
//         PROVIDES_NORMAL_DERIVATIVE
//     };

//     static const int NbDerivatives   = ((Type & FitScaleDer) ? 1 : 0 ) + ((Type & FitSpaceDer) ? DataPoint::Dim : 0);
//     static const int DerStorageOrder = (Type & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor;

// public:
//     typedef typename Base::Scalar     Scalar;     /*!< \brief Inherited scalar type*/
//     typedef typename Base::VectorType VectorType; /*!< \brief Inherited vector type*/
//     typedef typename Base::MatrixType MatrixType; /*!< \brief Inherited matrix type*/
//     typedef typename Base::WFunctor   WFunctor;   /*!< \brief Weight Function*/

//     /*! \brief Static array of scalars with a size adapted to the differentiation type */
//     typedef Eigen::Matrix<Scalar, DataPoint::Dim, NbDerivatives, DerStorageOrder> VectorArray;

//     /*! \brief Static array of scalars with a size adapted to the differentiation type */
//     typedef Eigen::Matrix<Scalar, 1, NbDerivatives> ScalarArray;
// private:
//     // computation data
//     ScalarArray m_dSumW;      /*!< \brief Sum of weight derivatives */
//     MatrixType  m_dCov[NbDerivatives];

//     VectorArray m_dCog;       /*!< \brief Derivatives of the centroid */
//     VectorArray m_dNormal;    /*!< \brief Derivatives of the hyper-surface normal */
//     ScalarArray m_dDist;      /*!< \brief Derivatives of the MLS scalar field */

// public:

//     /************************************************************************/
//     /* Initialization                                                       */
//     /************************************************************************/
//     /*! \see Concept::FittingProcedureConcept::init() */
//     PONCA_MULTIARCH void init(const VectorType &evalPos);

//     /************************************************************************/
//     /* Processing                                                           */
//     /************************************************************************/
//     /*! \see Concept::FittingProcedureConcept::addNeighbor() */
//     PONCA_MULTIARCH bool addNeighbor(const DataPoint  &nei);
//     /*! \see Concept::FittingProcedureConcept::finalize() */
//     PONCA_MULTIARCH FIT_RESULT finalize();


//     /**************************************************************************/
//     /* Use results                                                            */
//     /**************************************************************************/

//     /*! \brief Returns the derivatives of the scalar field at the evaluation point */
//     PONCA_MULTIARCH inline ScalarArray dPotential() const { return m_dDist; }

//     /*! \brief Returns the derivatives of the primitive normal */
//     PONCA_MULTIARCH inline VectorArray dNormal() const { return m_dNormal; }

//     /*! \brief State specified at compilation time to differenciate the fit in scale */
//     PONCA_MULTIARCH inline bool isScaleDer() const {return bool(Type & FitScaleDer);}
//     /*! \brief State specified at compilation time to differenciate the fit in space */
//     PONCA_MULTIARCH inline bool isSpaceDer() const {return bool(Type & FitSpaceDer);}
//     /*! \brief Number of dimensions used for the differentiation */
//     PONCA_MULTIARCH inline unsigned int derDimension() const { return NbDerivatives;}

// }; //class CovarianceSurfaceDer

// }// namespace internal

// /*!
//     \brief Differentiation in scale of the LeastSquareSurfaceFit
//     \inherit Concept::FittingExtensionConcept

//     Requierement:
//     \verbatim PROVIDES_COVARIANCE_SURFACE \endverbatim
//     Provide:
//     \verbatim PROVIDES_COVARIANCE_SURFACE_SCALE_DERIVATIVE \endverbatim

//     \ingroup fitting
// */
// template < class DataPoint, class _WFunctor, typename T>
// class CovarianceSurfaceScaleDer:public internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitScaleDer>
// {
// protected:
//     /*! \brief Inherited class */
//     typedef internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitScaleDer> Base;
//     enum { PROVIDES_COVARIANCE_SURFACE_SCALE_DERIVATIVE, PROVIDES_NORMAL_SCALE_DERIVATIVE };
// };


// /*!
//     \brief Spatial differentiation of the LeastSquareSurfaceFit
//     \inherit Concept::FittingExtensionConcept

//     Requierement:
//     \verbatim PROVIDES_COVARIANCE_SURFACE \endverbatim
//     Provide:
//     \verbatim PROVIDES_COVARIANCE_SURFACE_SPACE_DERIVATIVE \endverbatim

//     \ingroup fitting
// */
// template < class DataPoint, class _WFunctor, typename T>
// class CovarianceSurfaceSpaceDer:public internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitSpaceDer>
// {
// protected:
//     /*! \brief Inherited class */
//     typedef internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitSpaceDer> Base;
//     enum { PROVIDES_COVARIANCE_SURFACE_SPACE_DERIVATIVE, PROVIDES_NORMAL_SPACE_DERIVATIVE };
// };


// /*!
//     \brief Differentiation both in scale and space of the LeastSquareSurfaceFit
//     \inherit Concept::FittingExtensionConcept

//     Requierement:
//     \verbatim PROVIDES_COVARIANCE_SURFACE \endverbatim
//     Provide:
//     \verbatim PROVIDES_COVARIANCE_SURFACE_SCALE_DERIVATIVE
//     PROVIDES_COVARIANCE_SURFACE_SPACE_DERIVATIVE
//     \endverbatim

//     \ingroup fitting
// */
// template < class DataPoint, class _WFunctor, typename T>
// class CovarianceSurfaceScaleSpaceDer:public internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer>
// {
// protected:
//     /*! \brief Inherited class */
//     typedef internal::CovarianceSurfaceDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer> Base;
//     enum
//     {
//         PROVIDES_COVARIANCE_SURFACE_SCALE_DERIVATIVE,
//         PROVIDES_COVARIANCE_SURFACE_SPACE_DERIVATIVE,
//         PROVIDES_NORMAL_SCALE_DERIVATIVE,
//         PROVIDES_NORMAL_SPACE_DERIVATIVE
//     };
// };

#include "leastSquareSurfaceFit.hpp"

} //namespace Ponca
