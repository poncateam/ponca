/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./primitive.h"

namespace Ponca
{

/*!
    \brief Compute the barycenter of the input points
    \inherit Concept::FittingProcedureConcept

    \warning The barycenter is not stored explicitly, but rather computed from the sum of the neighbors positions.

    \todo Add scale and space derivatives

    This primitive provides:
    \verbatim PROVIDES_MEAN_POSITION \endverbatim

    \ingroup fitting
*/
    template < class DataPoint, class _WFunctor, typename T>
    class MeanPosition : public T
    {
    private:
        using Base = T;

    protected:
        enum
        {
            PROVIDES_MEAN_POSITION
        };

    public:
        using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
        using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
        using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/

    protected:
        VectorType m_sumP; /*!< \brief Sum of the input points vectors */
        Scalar     m_sumW; /*!< \brief Sum of queries weight */
        WFunctor   m_w;    /*!< \brief Weight function (must inherits BaseWeightFunc) */

    public:
        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MeanPosition() :
            Base(), m_sumP( VectorType::Zero() ), m_sumW(Scalar(1)) {}

        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
        PONCA_MULTIARCH inline void setWeightFunc (const WFunctor& _w) { m_w  = _w; }

        /*! \copydoc Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init (const VectorType& _evalPos);

        /**************************************************************************/
        /* Processing                                                             */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
        PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei);

        /*! \copydoc Concept::FittingProcedureConcept::finalize() */
        PONCA_MULTIARCH inline FIT_RESULT finalize();

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
        //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
        PONCA_MULTIARCH inline VectorType barycenter () const
        {
            return (m_sumP / m_sumW);
        }

    }; //class MeanPosition

/*!
    \brief Compute the barycenter of the input points + their normals
    \inherit Concept::FittingProcedureConcept

    \warning The barycenter is not stored explicitly, but rather computed from the sum of the neighbors positions and
    normals.

    This primitive provides:
    \verbatim PROVIDES_MEAN_NORMAL \endverbatim

    \see MeanPosition

    \note This class should not derive from MeanPosition, as we might want to compute mean normals but without mean
    positions. This is done this way currently, because we do not want to duplicate the weighting functor, which is
    currently stored in MeanPosition.

    \todo Add scale and space derivatives

    \ingroup fitting
*/
    template < class DataPoint, class _WFunctor, typename T >
    class MeanNormal : public MeanPosition<DataPoint, _WFunctor, T>
    {
    private:
        using Base = MeanPosition<DataPoint, _WFunctor, T>;

    protected:
        enum
        {
            PROVIDES_MEAN_NORMAL
        };

    public:
        using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
        using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
        using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/

    protected:
        VectorType  m_sumN;    /*!< \brief Sum of the normal vectors */

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MeanNormal() : Base() {}

        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/

        /*! \copydoc Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init (const VectorType& _evalPos);

        /**************************************************************************/
        /* Processing                                                             */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
        PONCA_MULTIARCH inline bool addNeighbor(const DataPoint &_nei);
    }; //class BarycenterWithNormal

#include "mean.hpp"

} //namespace Ponca
