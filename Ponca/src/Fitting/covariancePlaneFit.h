/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./plane.h"
#include "./primitive.h"

#include <Eigen/Eigenvalues>

namespace Ponca
{

/*!
    \brief Plane fitting procedure using only points position

    This class can also computes the surface variation measure introduced in
    \cite Pauly:2002:PSSimplification. The solver used to analyse the covariance
    matrix is stored for further use.

    \inherit Concept::FittingProcedureConcept

    \see Plane
    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T >
class CovariancePlaneFit : public MeanPosition<DataPoint, _WFunctor,
                                  Plane<DataPoint, _WFunctor>>
{
private:
    using Base = MeanPosition<DataPoint, _WFunctor, Plane<DataPoint, _WFunctor>>;

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE &&
                Base::PROVIDES_MEAN_POSITION,
        /*!
         * \brief Expose a method worldToTangentPlane(VectorType), which turns a point
         * in ambient 3D space to the tangent plane.
         * \see worldToTangentPlane
         * \see tangentPlaneToWorld
         */
        PROVIDES_TANGENT_PLANE_BASIS
    };

public:
    using Scalar     = typename Base::Scalar;          /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType;      /*!< \brief Inherited vector type*/
    using MatrixType = typename DataPoint::MatrixType; /*!< \brief Inherited matrix type*/
    using WFunctor   = typename Base::WFunctor;        /*!< \brief Weight Function*/
    /*! \brief Solver used to analyse the covariance matrix*/
    using Solver     = Eigen::SelfAdjointEigenSolver<MatrixType>;


 protected:
    // computation data
    MatrixType m_cov {MatrixType::Zero()};     /*!< \brief Covariance matrix */

    Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline CovariancePlaneFit() = default;

    /*! \brief Explicit conversion to CovariancePlaneFit, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    CovariancePlaneFit<DataPoint, WFunctor, T>& covariancePlaneFit()
    { return * static_cast<CovariancePlaneFit<DataPoint, WFunctor, T>*>(this); }

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

    /**************************************************************************/
    /* Results                                                                */
    /**************************************************************************/

    using Base::potential;

    /*! \brief Reading access to the Solver used to analyse the covariance
      matrix */
    PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }

    /*! \brief Implements \cite Pauly:2002:PSSimplification surface variation.

        It computes the ratio \f$ d \frac{\lambda_0}{\sum_i \lambda_i} \f$ with \c d the dimension of the ambient space.

        \return 0 for invalid fits
    */
    PONCA_MULTIARCH inline Scalar surfaceVariation() const;

    /*!
     * \brief Express a point in ambient space relatively to the tangent plane.
     * Output vector is: [h, u, v]^T, where u, v are 2d coordinates on the plane,
     * and h the height of the sample.
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType worldToTangentPlane(const VectorType &_q) const;

    /*!
     * \brief Transform a point from the tangent plane [h, u, v]^T to ambient space
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType tangentPlaneToWorld(const VectorType &_q) const;
}; //class CovariancePlaneFit

namespace internal {

using ::Ponca::internal::FitSpaceDer;
using ::Ponca::internal::FitScaleDer;

/*!
    \brief Internal generic class computing the derivatives of covariance plane fits
    \inherit Concept::FittingExtensionConcept

    The differentiation can be done automatically in scale and/or space, by
    combining the enum values FitScaleDer and FitSpaceDer in the template
    parameter Type.

    The differenciated values are stored in static arrays. The size of the
    arrays is computed with respect to the derivation type (scale and/or space)
    and the number of the dimension of the ambiant space.
    By convention, the scale derivatives are stored at index 0 when Type
    contains at least FitScaleDer. The size of these arrays can be known using
    derDimension(), and the differentiation type by isScaleDer() and
    isSpaceDer().

    \ingroup fitting
    \warning This class cannot be used directly, see CovariancePlaneScaleDer, CovariancePlaneSpaceDer and CovariancePlaneScaleSpaceDer
*/
template < class DataPoint, class _WFunctor, typename T, int Type>
class CovariancePlaneDer : public MeanPositionDer<DataPoint, _WFunctor, T, Type>
{
private:
    using Base = MeanPositionDer<DataPoint, _WFunctor, T, Type>; /*!< \brief Generic base type */


protected:
    enum
    {
        Check = Base::PROVIDES_PLANE,             /*!< \brief Needs plane */
        PROVIDES_COVARIANCE_PLANE_DERIVATIVE,      /*!< \brief Provides derivatives for hyper-planes */
        PROVIDES_NORMAL_DERIVATIVE
    };


public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using MatrixType = typename Base::MatrixType; /*!< \brief Inherited matrix type*/
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/
    using ScalarArray = typename Base::ScalarArray;
    using VectorArray = typename Base::VectorArray;

private:
    // computation data
    MatrixType  m_dCov[Base::NbDerivatives];

    VectorArray m_dNormal;    /*!< \brief Derivatives of the hyper-plane normal */
    ScalarArray m_dDist;      /*!< \brief Derivatives of the MLS scalar field */

public:

    /************************************************************************/
    /* Initialization                                                       */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::init() */
    PONCA_MULTIARCH void init(const VectorType &evalPos);

    /************************************************************************/
    /* Processing                                                           */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::addLocalNeighbor() */
    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);
    /*! \see Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH FIT_RESULT finalize();


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/

    /*! \brief Returns the derivatives of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline ScalarArray dPotential() const { return m_dDist; }

    /*! \brief Returns the derivatives of the primitive normal */
    PONCA_MULTIARCH inline VectorArray dNormal() const { return m_dNormal; }

}; //class CovariancePlaneDer

}// namespace internal

/*!
    \brief Differentiation in scale of the CovariancePlaneFit
    \inherit Concept::FittingExtensionConcept

    Requirements:
    \verbatim PROVIDES_COVARIANCE_PLANE \endverbatim
    Provides:
    \verbatim PROVIDES_COVARIANCE_PLANE_SCALE_DERIVATIVE \endverbatim

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T>
class CovariancePlaneScaleDer:public internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitScaleDer>
{
protected:
    /*! \brief Inherited class */
    typedef internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitScaleDer> Base;
    enum { PROVIDES_COVARIANCE_PLANE_SCALE_DERIVATIVE, PROVIDES_NORMAL_SCALE_DERIVATIVE };
};


/*!
    \brief Spatial differentiation of the CovariancePlaneFit
    \inherit Concept::FittingExtensionConcept

    Requirements:
    \verbatim PROVIDES_COVARIANCE_PLANE \endverbatim
    Provides:
    \verbatim PROVIDES_COVARIANCE_PLANE_SPACE_DERIVATIVE \endverbatim

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T>
class CovariancePlaneSpaceDer:public internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitSpaceDer>
{
protected:
    /*! \brief Inherited class */
    typedef internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitSpaceDer> Base;
    enum { PROVIDES_COVARIANCE_PLANE_SPACE_DERIVATIVE, PROVIDES_NORMAL_SPACE_DERIVATIVE };
};


/*!
    \brief Differentiation both in scale and space of the CovariancePlaneFit
    \inherit Concept::FittingExtensionConcept

    Requirements:
    \verbatim PROVIDES_COVARIANCE_PLANE \endverbatim
    Provides:
    \verbatim PROVIDES_COVARIANCE_PLANE_SCALE_DERIVATIVE
    PROVIDES_COVARIANCE_PLANE_SPACE_DERIVATIVE
    \endverbatim

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T>
class CovariancePlaneScaleSpaceDer:public internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer>
{
protected:
    /*! \brief Inherited class */
    typedef internal::CovariancePlaneDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer> Base;
    enum
    {
        PROVIDES_COVARIANCE_PLANE_SCALE_DERIVATIVE,
        PROVIDES_COVARIANCE_PLANE_SPACE_DERIVATIVE,
        PROVIDES_NORMAL_SCALE_DERIVATIVE,
        PROVIDES_NORMAL_SPACE_DERIVATIVE
    };
};

#include "covariancePlaneFit.hpp"

} //namespace Ponca
