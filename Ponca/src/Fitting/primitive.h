/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "enums.h"

namespace Ponca
{

/*!
    \brief Primitive base class.

    This class stores and provides public access to the fitting state, and must
    be inherited by classes implementing new primitives.

    Protected fields #m_eCurrentState and #m_nbNeighbors should be updated
    during the fitting process by the inheriting class.

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T = void  >
class PrimitiveBase
{
public:
    using Scalar     = typename DataPoint::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename DataPoint::VectorType; /*!< \brief Inherited vector type*/
    using WFunctor   = _WFunctor;                      /*!< \brief Weight Function*/

protected:

    //! \brief Represent the current state of the fit (finalize function
    //! update the state)
    FIT_RESULT m_eCurrentState {UNDEFINED};

    //! \brief Give the number of neighbors
    int m_nbNeighbors {0};

    //! \brief Weight function (must inherits BaseWeightFunc)
    WFunctor   m_w;

    //! \brief Sum of the neighbors weights
    Scalar m_sumW {0};

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline PrimitiveBase() = default;

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
    PONCA_MULTIARCH inline void setWeightFunc (const WFunctor& _w) { m_w  = _w; }

    /*! \brief Reset fitting state status */
    PONCA_MULTIARCH inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        m_eCurrentState = UNDEFINED;
        m_nbNeighbors = 0;
        m_sumW = Scalar(0);
        m_w.init( _basisCenter );
    }

    /*! \brief Is the primitive well fitted an ready to use (finalize has been
    called)
    \warning The fit can be unstable (having neighbors between 3 and 6) */
    PONCA_MULTIARCH inline bool isReady() const
    {
        return (m_eCurrentState == STABLE) || (m_eCurrentState == UNSTABLE);
    }

    /*! \brief Is the plane fitted an ready to use (finalize has been called
    and the result is stable, eq. having more than 6 neighbors) */
    PONCA_MULTIARCH inline bool isStable() const { return m_eCurrentState == STABLE; }

    /*! \return the current test of the fit */
    PONCA_MULTIARCH inline FIT_RESULT getCurrentState() const
    {
        return m_eCurrentState;
    }

    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar, const VectorType &, const DataPoint &) {
        return true;
    }

    PONCA_MULTIARCH inline FIT_RESULT finalize(){
        // handle specific configurations
        // We need to have at least one neighbor to compute the mean
        if (m_sumW == Scalar(0.) || m_nbNeighbors < 1) {
            init( m_w.basisCenter() );
            return m_eCurrentState = UNDEFINED;
        }
        return m_eCurrentState = STABLE;
    }

}; //class Plane

}
