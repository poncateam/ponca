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

        /// \brief Barycenter of the input points
        ///
        /// Defined as \f$ b(\mathbf{x}) = \frac{\sum_i w_\mathbf{x}(\mathbf{p_i}) \mathbf{p_i}}{\sum_i w_\mathbf{x}(\mathbf{p_i})} \f$,
        ///  where \f$\left[\mathbf{p_i} \in \text{neighborhood}(\mathbf{x})\right]\f$ are all the point samples in \f$\mathbf{x}\f$'s neighborhood
        PONCA_MULTIARCH inline VectorType barycenter() const {
            return (m_sumP / Base::getWeightSum());
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

        /// \brief Compute derivatives of the barycenter.
        /// \see MeanPosition::barycenter()
        ///
        /// ### Step-by-step derivation from the barycenter definition
        /// Given the definition of the barycenter \f$ b(\mathbf{x}) = \frac{\sum_i w_\mathbf{x}(\mathbf{p_i}) \mathbf{p_i}}{\sum_i w_\mathbf{x}(\mathbf{p_i})} \f$,
        /// where \f$\left[\mathbf{p_i} \in \text{neighborhood}(\mathbf{x})\right]\f$ are all the point samples in \f$\mathbf{x}\f$'s neighborhood.
        ///
        /// We denote \f$ t(\mathbf{x}) = \sum_i w_\mathbf{x}(\mathbf{p_i}) \mathbf{p_i} \f$ and \f$ s(\mathbf{x}) = \sum_i w_\mathbf{x}(\mathbf{p_i})\f$,
        /// such that \f$ b(\mathbf{x}) = \frac{t(\mathbf{x})}{s(\mathbf{x})}\f$.
        ///
        /// By definition, \f$ b'(\mathbf{x}) = \frac{s(\mathbf{x})t'(\mathbf{x}) - t(\mathbf{x})s'(\mathbf{x})}{s(\mathbf{x})^2}\f$.
        /// We have \f$ s'(\mathbf{x}) = \sum_i w'_\mathbf{x}(\mathbf{p_i}) \f$.
        ///
        /// We rewrite \f$ t(\mathbf{x}) = \sum u(\mathbf{x})v(\mathbf{x}) \f$, with \f$ u(\mathbf{x}) = w_\mathbf{x}(\mathbf{p_i}) \f$ and \f$ v(\mathbf{x}) = \mathbf{p_i} \f$.
        ///
        /// As the point cloud coordinates are constants, \f$v(\mathbf{x})\f$ is constant, its derivative is null, and so \f$  t'(\mathbf{x}) = \sum_i u'(\mathbf{x}) v(\mathbf{x}) = \sum_i w'_\mathbf{x}(\mathbf{p_i}) \mathbf{p_i} \f$.
        ///
        /// Which leads to \f$ b'(\mathbf{x}) = \frac{\sum_i w'_\mathbf{x}(\mathbf{p_i}) \mathbf{p_i} - b(\mathbf{x})\sum w'(\mathbf{x})}{\sum_i w_\mathbf{x}(\mathbf{p_i})} \f$
        ///
        /// \note This code is not directly tested, but rather indirectly by testing CovariancePlaneDer::dNormal()
        PONCA_MULTIARCH VectorArray barycenterDerivatives() const
        {
            return ( m_dSumP - Base::barycenter() * Base::m_dSumW ) / Base::getWeightSum();
        }

    }; //class MeanPositionDer

#include "mean.hpp"

} //namespace Ponca
