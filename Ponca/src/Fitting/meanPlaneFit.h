/*
 Copyright (C) 2018 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./plane.h"
#include "./mean.h"
#include "./localFrame.h"

namespace Ponca
{

/*!
    \brief Plane fitting procedure computing the frame plane using mean position and normal

    \inherit Concept::FittingProcedureConcept

    \see Plane
    \see PlaneFrame
*/
template < class DataPoint, class _NFilter, typename T >
class MeanPlaneFitImpl : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME && Base::PROVIDES_PLANE };


public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlaneFitImpl,meanPlaneFit)

    /*!
     * \brief Fitting of the plane using the mean position and mean normal
     * We use the localFrame class to sorte the frame information.
     * Given the mean normal, we can compute the frame plane.
     * m_u and m_v are computed using the cross product, to ensure orthogonality.
     * \see LocalFrame
     */
    PONCA_FITTING_APIDOC_FINALIZE
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        // handle specific configurations
        if(Base::finalize() == STABLE)
        {
            if (Base::plane().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
            Base::setPlane(Base::m_sumN / Base::m_sumW, Base::barycenter());
            VectorType norm = Base::plane().normal();
            VectorType a;
            if (std::abs(norm.x()) > std::abs(norm.z())) {
                a = VectorType(-norm.y(), norm.x(), 0);
            } else {
                a = VectorType(0, -norm.z(), norm.y());
            }
            a.normalize();
            // use cross product to generate a orthogonal basis
            VectorType m_u = norm.cross(a);
            m_u.normalize();
            VectorType m_v = norm.cross(m_u);
            m_v.normalize();
            Base::setFrameUV (m_u, m_v);
        }
        return Base::m_eCurrentState;
    }
    PONCA_FITTING_IS_SIGNED(true)
}; //class MeanPlaneFitImpl

/// \brief Helper alias for Plane fitting on points using MeanPlaneFitImpl
//! [MeanPlaneFit Definition]
    template < class DataPoint, class _NFilter, typename T>
    using MeanPlaneFit =
    MeanPlaneFitImpl<DataPoint, _NFilter,
        MeanNormal<DataPoint, _NFilter,
            MeanPosition<DataPoint, _NFilter,
                LocalFrame<DataPoint, _NFilter,
                    Plane<DataPoint, _NFilter, T>>>>>;
//! [MeanPlaneFit Definition]

} //namespace Ponca
