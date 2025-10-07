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

    \note This class should not be inherited explicitly: this is done by the
    #Basket class.
*/
template < class DataPoint, class _NFilter, typename T = void  >
class PrimitiveBase
{
protected:
    enum {
        PROVIDES_PRIMITIVE_BASE,    /*!< \brief Provides base API for primitives*/
    };

public:
    using Scalar         = typename DataPoint::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType     = typename DataPoint::VectorType; /*!< \brief Inherited vector type*/
    using NeighborFilter = _NFilter;                       /*!< \brief Filter applied on each neighbor*/

private:
    //! \brief Number of neighbors
    int m_nbNeighbors {0};

    //! \brief Sum of the neighbors weights
    Scalar m_sumW {0};

    //! \brief Neighborhood filter
    NeighborFilter   m_nFilter;

protected:

    //! \brief Represent the current state of the fit (finalize function
    //! update the state)
    FIT_RESULT m_eCurrentState {UNDEFINED};

public:
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    PONCA_FITTING_APIDOC_SETWFUNC
    PONCA_MULTIARCH inline void setNeighborFilter (const NeighborFilter& _nFilter) {
        m_nFilter  = _nFilter;
    }

    PONCA_FITTING_APIDOC_INIT
    PONCA_MULTIARCH inline void init()
    {
        m_eCurrentState = UNDEFINED;
        startNewPass();
    }

    /*! \brief Is the primitive well fitted and ready to use (finalize has been
    called) ?
    \warning The fit can be unstable (having neighbors between 3 and 6) */
    PONCA_MULTIARCH inline bool isReady() const
    {
        return (m_eCurrentState == STABLE) || (m_eCurrentState == UNSTABLE);
    }

    /*! \brief Is the fitted primitive ready to use (finalize has been called and the result is stable) */
    PONCA_MULTIARCH inline bool isStable() const { return m_eCurrentState == STABLE; }

    /*! \brief Is another pass required for fitting (finalize has been called and the result is #NEED_OTHER_PASS)
     * \see startNewPass */
    PONCA_MULTIARCH inline bool needAnotherPass() const { return m_eCurrentState == NEED_OTHER_PASS; }

    /*! \brief Get number of points added in the neighborhood (with non negative weight)  */
    PONCA_MULTIARCH inline int getNumNeighbors() const { return m_nbNeighbors; }

    /*! \brief Get the sum of the weights */
    PONCA_MULTIARCH inline Scalar getWeightSum() const { return m_sumW; }

    /*! \brief To be called when starting a new processing pass, ie. when `getCurrentState()==#NEED_ANOTHER_PASS` */
    PONCA_MULTIARCH inline void startNewPass() {
        m_nbNeighbors = 0;
        m_sumW = Scalar(0);
    }

    /*! \brief Read access to the NeighborFilter \see setNeighborFilter */
    PONCA_MULTIARCH inline const NeighborFilter& getNeighborFilter() const
    {
        return m_nFilter;
    }


    /*! \return the current test of the fit */
    PONCA_MULTIARCH inline FIT_RESULT getCurrentState() const
    {
        return m_eCurrentState;
    }

    PONCA_FITTING_APIDOC_ADDNEIGHBOR
    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &, const DataPoint &) {
        m_sumW += w;
        ++(m_nbNeighbors);
        return true;
    }

    PONCA_FITTING_APIDOC_FINALIZE
    PONCA_MULTIARCH inline FIT_RESULT finalize(){
        // handle specific configurations
        if (m_sumW == Scalar(0.) || m_nbNeighbors < 1) {
            return m_eCurrentState = UNDEFINED;
        }
        return m_eCurrentState = STABLE;
    }

}; //class Primitive





/**
    \brief Generic class performing the Fit derivation
    \inherit Concept::FittingExtensionConcept

    The differentiation can be done automatically in scale and/or space, by
    combining the enum values FitScaleDer and FitSpaceDer in the template
    parameter Type.

    The differentiated values are stored in static arrays. The size of the
    arrays is computed with respect to the derivation type (scale and/or space)
    and the number of the dimension of the ambiant space.
    By convention, the scale derivatives are stored at index 0 when Type
    contains at least FitScaleDer. The size of these arrays can be known using
    derDimension(), and the differentiation type by isScaleDer() and
    isSpaceDer().

    Thanks to the BasketDiff definition, we know that PrimitiveDer has Primitive
    as base class (through the Basket). As a result, this class first asks to
    compute the Fit, and if it works properly, compute the weight derivatives.
 */
template < class DataPoint, class _NFilter, int Type, typename T>
class PrimitiveDer : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE
protected:
    enum {
        Check = Base::PROVIDES_PRIMITIVE_BASE,    /*!< \brief Provides base API for primitives*/
        PROVIDES_PRIMITIVE_DERIVATIVE
    };

protected:
    static constexpr  int NbDerivatives   = ((Type & FitScaleDer) ? 1 : 0 ) + ((Type & FitSpaceDer) ? DataPoint::Dim : 0);
    static constexpr  int DerStorageOrder = (Type & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor;

public:

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
    PONCA_MULTIARCH inline void init()
    { Base::init(); m_dSumW.setZero(); }

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
        if( Base::addLocalNeighbor(w, localQ, attributes) ) { // call the Primitive Fit (without dw)
            int spaceId = (Type & FitScaleDer) ? 1 : 0;
            // compute weight
            if (Type & FitScaleDer)
                dw[0] = Base::getNeighborFilter().scaledw(attributes.pos(), attributes);

            if (Type & FitSpaceDer)
                dw.template segment<int(DataPoint::Dim)>(spaceId) = -Base::getNeighborFilter().spacedw(attributes.pos(), attributes).transpose();

            m_dSumW += dw;
            return true;
        }
        return false;
    }

    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief State specified at compilation time to differenciate the fit in scale */
    PONCA_MULTIARCH inline constexpr bool isScaleDer() const {return bool(Type & FitScaleDer);}
    /*! \brief State specified at compilation time to differenciate the fit in space */
    PONCA_MULTIARCH inline constexpr bool isSpaceDer() const {return bool(Type & FitSpaceDer);}
    /*! \brief Number of dimensions used for the differentiation */
    PONCA_MULTIARCH inline constexpr unsigned int derDimension() const { return NbDerivatives;}

};



}
