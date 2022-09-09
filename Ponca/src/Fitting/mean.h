/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./primitive.h"

namespace Ponca {

/*!
    \brief Compute the barycenter of the input points
    \inherit Concept::FittingProcedureConcept

    \warning The barycenter is not stored explicitly, but rather computed from the sum of the neighbors positions.

    This primitive provides:
    \verbatim PROVIDES_MEAN_POSITION \endverbatim

    \ingroup fitting
*/
    template<class DataPoint, class _WFunctor, typename T>
    class MeanPosition : public T {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

    protected:
        enum { PROVIDES_MEAN_POSITION };
        VectorType m_sumP {VectorType::Zero()}; /*!< \brief Sum of the input points vectors */

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(MeanPosition,meanPosition)
        PONCA_FITTING_DECLARE_INIT
        PONCA_FITTING_DECLARE_ADDNEIGHBOR

        //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
        PONCA_MULTIARCH inline VectorType barycenter() const {
            return (m_sumP / Base::m_sumW);
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
    template<class DataPoint, class _WFunctor, typename T>
    class MeanNormal : public T {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

    protected:
        enum { PROVIDES_MEAN_NORMAL };
        VectorType m_sumN;    /*!< \brief Sum of the normal vectors */

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(MeanNormal,meanNormal)
        PONCA_FITTING_DECLARE_INIT
        PONCA_FITTING_DECLARE_ADDNEIGHBOR
    }; //class MeanNormal

    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    class MeanPositionDer : public T {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

    protected:
        enum {
            Check = Base::PROVIDES_PRIMITIVE_DERIVATIVE &&
                    Base::PROVIDES_MEAN_POSITION,
            PROVIDES_MEAN_POSITION_DERIVATIVE,    /*!< \brief Provides derivative of the mean position*/
        };

        /*! \brief Derivatives of the of the input points vectors */
        VectorArray m_dSumP {VectorArray::Zero()};

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(MeanPositionDer,meanPositionDer)
        PONCA_FITTING_DECLARE_INIT
        PONCA_FITTING_DECLARE_ADDNEIGHBOR_DER

        /// \brief Compute derivatives of the barycenter. \see MeanPosition::barycenter()
        PONCA_MULTIARCH VectorArray barycenterDerivatives() const
        {
            VectorArray barycenterDer;
            for(int k=0; k<Base::NbDerivatives; ++k)
                barycenterDer.col(k) = (m_dSumP.col(k) - Base::m_dSumW(k) * Base::m_sumP) / Base::m_sumW;
            return barycenterDer;
        }

    }; //class MeanPositionDer

#include "mean.hpp"

} //namespace Ponca
