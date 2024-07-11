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
*/
template < class DataPoint, class _WFunctor, typename T >
class MeanPlaneFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_MEAN_POSITION && Base::PROVIDES_MEAN_NORMAL && Base::PROVIDES_PLANE };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlaneFitImpl,meanPlaneFit)

    PONCA_FITTING_APIDOC_FINALIZE
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        // handle specific configurations
        if(Base::finalize() == STABLE)
        {
            if (Base::plane().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
            Base::setPlane(Base::m_sumN / Base::getWeightSum(), Base::barycenterLocal());
        }
        return Base::m_eCurrentState;
    }
    PONCA_FITTING_IS_SIGNED(true)
}; //class MeanPlaneFitImpl

/// \brief Helper alias for Plane fitting on points using MeanPlaneFitImpl
//! [MeanPlaneFit Definition]
    template < class DataPoint, class _WFunctor, typename T>
    using MeanPlaneFit =
    MeanPlaneFitImpl<DataPoint, _WFunctor,
        MeanNormal<DataPoint, _WFunctor,
            MeanPosition<DataPoint, _WFunctor,
                Plane<DataPoint, _WFunctor,T>>>>;
//! [MeanPlaneFit Definition]
} //namespace Ponca
