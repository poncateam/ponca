/*
 Copyright (C) 2018 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./plane.h"
#include "./mean.h"

namespace Ponca
{

/*!
    \brief Plane fitting procedure computing the mean position and orientation
    from oriented points

    \inherit Concept::FittingProcedureConcept

    \see Plane

    \todo Add local frame computation to enable PROVIDES_TANGENT_PLANE_BASIS

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T >
class MeanPlaneFitImpl : public T
{
private:
    using Base = T;

protected:
    enum
    {
        Check = Base::PROVIDES_MEAN_POSITION && Base::PROVIDES_MEAN_NORMAL && Base::PROVIDES_PLANE
    };

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using MatrixType = typename DataPoint::MatrixType;
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline MeanPlaneFitImpl() = default;

    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlaneFitImpl,meanPlaneFit)

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        // handle specific configurations
        if(Base::finalize() == STABLE)
        {
            if (Base::plane().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
            Base::setPlane(Base::m_sumN / Base::m_sumW, Base::barycenter());
        }
        return Base::m_eCurrentState;
    }
}; //class MeanPlaneFitImpl

/// \brief Helper alias for Plane fitting on points using MeanPlaneFitImpl
/// \ingroup fittingalias
    template < class DataPoint, class _WFunctor, typename T>
    using MeanPlaneFit =
    MeanPlaneFitImpl<DataPoint, _WFunctor,
            MeanNormal<DataPoint, _WFunctor,
            MeanPosition<DataPoint, _WFunctor,
            Plane<DataPoint, _WFunctor,T>>>>;
} //namespace Ponca
