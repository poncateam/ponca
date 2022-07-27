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
    private:
        using Base = T;

    protected:
        enum {
            PROVIDES_MEAN_POSITION
        };

    public:
        using Scalar = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
        using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
        using WFunctor = typename Base::WFunctor;   /*!< \brief Weight Function*/

    protected:
        VectorType m_sumP {VectorType::Zero()}; /*!< \brief Sum of the input points vectors */

    public:
        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MeanPosition() = default;

        PONCA_EXPLICIT_CAST_OPERATORS(MeanPosition,meanPosition)

        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/
        /*! \copydoc Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init(const VectorType &_evalPos);

        /**************************************************************************/
        /* Processing                                                             */
        /**************************************************************************/
        /*! \copydoc Concept::FittingExtensionConcept::addLocalNeighbor() */
        PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
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
    private:
        using Base = T;

    protected:
        enum {
            PROVIDES_MEAN_NORMAL
        };

    public:
        using Scalar = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
        using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
        using WFunctor = typename Base::WFunctor;   /*!< \brief Weight Function*/

    protected:
        VectorType m_sumN;    /*!< \brief Sum of the normal vectors */

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MeanNormal() : Base() {}

        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/

        /*! \copydoc Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init(const VectorType &_evalPos);

        /**************************************************************************/
        /* Processing                                                             */
        /**************************************************************************/
        /*! \copydoc Concept::FittingExtensionConcept::addLocalNeighbor() */
        PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);
    }; //class MeanNormal

    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    class MeanPositionDer : public T {
    private:
        using Base = T; /*!< \brief Generic base type */


    protected:
        enum {
            Check = Base::PROVIDES_PRIMITIVE_DERIVATIVE &&
                    Base::PROVIDES_MEAN_POSITION,
            PROVIDES_MEAN_POSITION_DERIVATIVE,    /*!< \brief Provides derivative of the mean position*/
        };

    public:
        using Scalar = typename Base::Scalar;
        using VectorType = typename Base::VectorType;
        using WFunctor = typename Base::WFunctor;
        using ScalarArray = typename Base::ScalarArray;
        using VectorArray = typename Base::VectorArray;

    protected:
        /*! \brief Derivatives of the of the input points vectors */
        VectorArray m_dSumP {VectorArray::Zero()};

    public:
        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MeanPositionDer() = default;

        PONCA_EXPLICIT_CAST_OPERATORS_DER(MeanPositionDer,meanPositionDer)
        /************************************************************************/
        /* Initialization                                                       */
        /************************************************************************/
        /*! \see Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH void init(const VectorType &evalPos);

        /************************************************************************/
        /* Processing                                                           */
        /************************************************************************/
        /*! \see Concept::FittingProcedureConcept::addLocalNeighbor() */
        PONCA_MULTIARCH inline bool
        addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes, ScalarArray &dw);

        /*! \see Concept::FittingProcedureConcept::finalize() */
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
