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
#include "./mean.h"          // used to define CovarianceLineFit
#include "./covarianceFit.h" // use to define CovariancePlaneFit

#include <Eigen/Eigenvalues>

namespace Ponca
{

/*!
    \brief Plane fitting procedure using only points position

    This class can also computes the surface variation measure introduced in
    \cite Pauly:2002:PSSimplification.

    \inherit Concept::FittingProcedureConcept

    \see Plane
    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T >
class CovariancePlaneFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE &&
                Base::PROVIDES_POSITION_COVARIANCE,
        /*!
         * \brief Expose a method worldToTangentPlane(VectorType), which turns a point
         * in ambient 3D space to the tangent plane.
         * \see worldToTangentPlane
         * \see tangentPlaneToWorld
         */
        PROVIDES_TANGENT_PLANE_BASIS
    };

public:
    /*! \brief Default constructor */
    PONCA_MULTIARCH inline CovariancePlaneFitImpl() = default;

    PONCA_EXPLICIT_CAST_OPERATORS(CovariancePlaneFitImpl,covariancePlaneFit)
    PONCA_FITTING_DECLARE_FINALIZE

    /**************************************************************************/
    /* Results                                                                */
    /**************************************************************************/

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
}; //class CovariancePlaneFitImpl

/// \brief Helper alias for Plane fitting on 3D points using CovariancePlaneFitImpl
/// \ingroup fittingalias
    template < class DataPoint, class _WFunctor, typename T>
    using CovariancePlaneFit =
    CovariancePlaneFitImpl<DataPoint, _WFunctor,
            CovarianceFitBase<DataPoint, _WFunctor,
                    MeanPosition<DataPoint, _WFunctor,
                            Plane<DataPoint, _WFunctor,T>>>>;

/*!
    \brief Internal generic class computing the derivatives of covariance plane fits
    \inherit Concept::FittingExtensionConcept

    \ingroup fitting
    \warning Defined in 3D only
*/
template < class DataPoint, class _WFunctor, int DiffType, typename T>
class CovariancePlaneDerImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE
    PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
    static_assert ( DataPoint::Dim == 3, "CovariancePlaneDer is only valid in 3D");

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE &
                Base::PROVIDES_POSITION_COVARIANCE_DERIVATIVE,
        PROVIDES_COVARIANCE_PLANE_DERIVATIVE,                    /*!< \brief Provides derivatives for hyper-planes */
        PROVIDES_NORMAL_DERIVATIVE
    };

private:
    VectorArray m_dNormal {VectorArray::Zero()};    /*!< \brief Derivatives of the hyper-plane normal */
    ScalarArray m_dDist {ScalarArray::Zero()};      /*!< \brief Derivatives of the MLS scalar field */

public:
    /*! \brief Default constructor */
    PONCA_MULTIARCH inline CovariancePlaneDerImpl() = default;

    PONCA_EXPLICIT_CAST_OPERATORS_DER(CovariancePlaneDerImpl,covariancePlaneDer)

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


/// \brief Helper alias for Plane fitting on 3D points using CovariancePlaneFitImpl
/// \ingroup fittingalias
template < class DataPoint, class _WFunctor, int DiffType, typename T>
using CovariancePlaneDer =
CovariancePlaneDerImpl<DataPoint, _WFunctor, DiffType,
        CovarianceFitDer<DataPoint, _WFunctor, DiffType,
                MeanPositionDer<DataPoint, _WFunctor, DiffType, T>>>;

#include "covariancePlaneFit.hpp"

} //namespace Ponca
