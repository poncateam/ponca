/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Eigenvalues>

namespace Ponca
{

/*!
    \brief Plane fitting procedure using only points position

    This class can also computes the surface variation measure introduced in
    \cite Pauly:2002:PSSimplification.

    \inherit Concept::FittingProcedureConcept

    \see Plane
*/
template < class DataPoint, class _NFilter, typename T >
class CovariancePlaneFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        /*!
         * \brief Fit the tangent plane and store it into Plane and LocalFrame which turn a point
         * in ambient 3D space to the tangent plane.
         * \see LocalFrame
         */
        Check = Base::PROVIDES_POSITION_COVARIANCE &&
                Base::PROVIDES_LOCAL_FRAME
        };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(CovariancePlaneFitImpl,covariancePlaneFit)
    PONCA_FITTING_DECLARE_FINALIZE
    PONCA_FITTING_IS_SIGNED(false)

}; //class CovariancePlaneFitImpl

/// \brief Helper alias for Plane fitting on 3D points using CovariancePlaneFitImpl
//! [CovariancePlaneFit Definition]
    template < class DataPoint, class _NFilter, typename T>
    using CovariancePlaneFit =
    CovariancePlaneFitImpl<DataPoint, _NFilter,
            CovarianceFitBase<DataPoint, _NFilter,
                    MeanPosition<DataPoint, _NFilter,
                        LocalFrame<DataPoint, _NFilter,
                            Plane<DataPoint, _NFilter,T>>>>>;
//! [CovariancePlaneFit Definition]

/*!
    \brief Internal generic class computing the derivatives of covariance plane fits
    \inherit Concept::FittingExtensionConcept

    \warning Defined in 3D only
*/
template < class DataPoint, class _NFilter, int DiffType, typename T>
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
    PONCA_EXPLICIT_CAST_OPERATORS_DER(CovariancePlaneDerImpl,covariancePlaneDer)

    /*! \see Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH FIT_RESULT finalize();

    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/

    // \brief Returns the derivatives of the scalar field at the evaluation point
    //! \see method `#isSigned` of the fit to check if the sign is reliable
    PONCA_MULTIARCH [[nodiscard]] inline ScalarArray dPotential() const { return m_dDist; }

    /*! \brief Returns the derivatives of the primitive normal */
    PONCA_MULTIARCH [[nodiscard]] inline VectorArray dNormal() const { return m_dNormal; }

}; //class CovariancePlaneDer


/// \brief Helper alias for Plane fitting on 3D points using CovariancePlaneFitImpl
template < class DataPoint, class _NFilter, int DiffType, typename T>
using CovariancePlaneDer =
CovariancePlaneDerImpl<DataPoint, _NFilter, DiffType,
        CovarianceFitDer<DataPoint, _NFilter, DiffType,
                MeanPositionDer<DataPoint, _NFilter, DiffType, T>>>;

#include "covariancePlaneFit.hpp"

} //namespace Ponca
