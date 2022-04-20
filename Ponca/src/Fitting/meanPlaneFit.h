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

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T >
class MeanPlaneFit : public MeanNormal<DataPoint, _WFunctor, Plane<DataPoint, _WFunctor>>
{
private:
    using Base = MeanNormal<DataPoint, _WFunctor, Plane<DataPoint, _WFunctor>>;

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE && Base::PROVIDES_MEAN_NORMAL
    };

public:

    /*! \brief Scalar type inherited from DataPoint*/
    typedef typename Base::Scalar     Scalar;
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename Base::VectorType VectorType;
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename Base::MatrixType MatrixType;

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline MeanPlaneFit() : Base() {}

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
        // handle specific configurations
        if(Base::finalize() == STABLE)
        {
            Base::setPlane(Base::m_sumN / Base::m_sumW, Base::barycenter());
        }

        return Base::m_eCurrentState;
    }
}; //class MeanPlaneFit
} //namespace Ponca
