/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include "./defines.h"
#include "./linePrimitive.h" // used to define CovarianceLineFit
#include "./mean.h"          // used to define CovarianceLineFit
#include "./covarianceFit.h" // used to define CovarianceLineFit

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
class CovarianceLineFitImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_LINE &&
                Base::PROVIDES_POSITION_COVARIANCE,
    };

public:
     /*! \brief Default constructor */
    PONCA_MULTIARCH inline CovarianceLineFitImpl() = default;

    PONCA_EXPLICIT_CAST_OPERATORS(CovarianceLineFitImpl,covarianceLineFit)

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        static const int smallestEigenValue = DataPoint::Dim - 1;
        if (Base::finalize() == STABLE)
            Base::setLine(Base::barycenter(), Base::m_solver.eigenvectors().col(smallestEigenValue).normalized());
        return Base::m_eCurrentState;
    }
};

/// \brief Helper alias for Line fitting on 3D points using CovarianceLineFitImpl
/// \ingroup fittingalias
template < class DataPoint, class _WFunctor, typename T>
using CovarianceLineFit =
                CovarianceLineFitImpl<DataPoint, _WFunctor,
                CovarianceFitBase<DataPoint, _WFunctor,
                MeanPosition<DataPoint, _WFunctor,
                Line<DataPoint, _WFunctor,T>>>>;

} //namespace Ponca
