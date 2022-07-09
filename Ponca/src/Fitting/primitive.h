/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "enums.h"
#include <Eigen/Dense>

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

}; //class Primitive


namespace internal
{
template < class DataPoint, class _WFunctor, typename T, int Type>
class PrimitiveDer : public T
{
private:
    typedef T Base; /*!< \brief Generic base type */

protected:
    static const int NbDerivatives   = ((Type & FitScaleDer) ? 1 : 0 ) + ((Type & FitSpaceDer) ? DataPoint::Dim : 0);
    static const int DerStorageOrder = (Type & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor;

public:
    using Scalar     = typename Base::Scalar;          /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType;      /*!< \brief Inherited vector type*/
    using MatrixType = typename DataPoint::MatrixType; /*!< \brief Inherited matrix type*/
    using WFunctor   = typename Base::WFunctor;        /*!< \brief Weight Function*/

    /*! \brief Static array of scalars with a size adapted to the differentiation type */
    typedef Eigen::Matrix<Scalar, DataPoint::Dim, NbDerivatives, DerStorageOrder> VectorArray;

    /*! \brief Static array of scalars with a size adapted to the differentiation type */
    typedef Eigen::Matrix<Scalar, 1, NbDerivatives> ScalarArray;

protected:
    // computation data
    ScalarArray m_dSumW;      /*!< \brief Sum of weight derivatives */

public:

    /************************************************************************/
    /* Initialization                                                       */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::init() */
    PONCA_MULTIARCH inline void init(const VectorType &_evalPos)
    { Base::init(_evalPos); m_dSumW.setZero(); }

    /************************************************************************/
    /* Processing                                                           */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::addLocalNeighbor()
     * \todo update doc (takes derivatives) */
    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w,
                                                 const VectorType &localQ,
                                                 const DataPoint &attributes,
                                                 ScalarArray &dw)
    {
        if( Base::addLocalNeighbor(w, localQ, attributes) ) {
            int spaceId = (Type & FitScaleDer) ? 1 : 0;
            // compute weight
            if (Type & FitScaleDer)
                dw[0] = Base::m_w.scaledw(attributes.pos(), attributes);

            if (Type & FitSpaceDer)
                dw.template segment<int(DataPoint::Dim)>(spaceId) = -Base::m_w.spacedw(attributes.pos(), attributes).transpose();

            m_dSumW += dw;
            return true;
        }
        return false;
    }

    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief State specified at compilation time to differenciate the fit in scale */
    PONCA_MULTIARCH inline bool isScaleDer() const {return bool(Type & FitScaleDer);}
    /*! \brief State specified at compilation time to differenciate the fit in space */
    PONCA_MULTIARCH inline bool isSpaceDer() const {return bool(Type & FitSpaceDer);}
    /*! \brief Number of dimensions used for the differentiation */
    PONCA_MULTIARCH inline unsigned int derDimension() const { return NbDerivatives;}

};
}
}
