/*
 Copyright (C) 2018 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./mean.h"       // used to define MeanPlaneFit
#include "./plane.h"      // used to define MeanPlaneFit
#include "./localFrame.h" // used to define MeanPlaneFit


namespace Ponca
{

/*!
    \brief Plane fitting procedure computing the frame plane using mean position and normal

    \inherit Concept::FittingProcedureConcept

    \see Plane
    \see localFrame
*/
template < class DataPoint, class _WFunctor, typename T >
class MeanPlaneFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_PLANE
                && Base::PROVIDES_MEAN_POSITION 
                && Base::PROVIDES_MEAN_NORMAL 
                && Base::PROVIDES_LOCAL_FRAME
         };


public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlaneFitImpl,meanPlaneFit)

    /*!
     * \brief This function fits the plane using mean normal and position.

     * We use the localFrame class to store the frame informations.
     * Given the mean normal, we can compute the frame plane.
     * m_u and m_v are computed using the cross product, to ensure orthogonality.
     * \see LocalFrame
     * \see computeFrameFromNormalVector
     */
    PONCA_FITTING_APIDOC_FINALIZE
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        // handle specific configurations
        if(Base::finalize() == STABLE)
        {
            if (Base::plane().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
            VectorType norm = Base::m_sumN / Base::getWeightSum();
            Base::setPlane(norm, Base::barycenter());
            Base::computeFrameFromNormalVector(norm);
        }
        return Base::m_eCurrentState;
    }

}; //class MeanPlaneFitImpl

/// \brief Helper alias for Plane fitting on points using MeanPlaneFitImpl
//! [MeanPlaneFit Definition]
    template < class DataPoint, class _WFunctor, typename T>
    using MeanPlaneFit =
    MeanPlaneFitImpl<DataPoint, _WFunctor,
        MeanNormal<DataPoint, _WFunctor,
            MeanPosition<DataPoint, _WFunctor,
                LocalFrame<DataPoint, _WFunctor,
                    Plane<DataPoint, _WFunctor,T>>>>>;
//! [MeanPlaneFit Definition]

} //namespace Ponca
