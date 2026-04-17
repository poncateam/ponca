/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "concepts.h"
#include "enums.h"
#include <Eigen/Dense>

namespace Ponca
{

    /*!
        \brief Base class of any computation unit of a Basket

        It is the first unit of the Basket class. It is used to define default fallback functions, and to store
        and provide public access to the neighborhood statistics, fitting status, ... It also stores the shared
        datastructures such as the NeighborFilter.
    */
    template <class DataPoint, class _NFilter, typename T = void>
    class BasketUnitBase
    {
    public:
        using Scalar         = typename DataPoint::Scalar;     /*!< \brief Inherited scalar type*/
        using VectorType     = typename DataPoint::VectorType; /*!< \brief Inherited vector type*/
        using NeighborFilter = _NFilter;                       /*!< \brief Filter applied on each neighbor*/

    private:
        //! \brief Number of neighbors
        int m_nbNeighbors{0};

        //! \brief Sum of the neighbors weights
        Scalar m_sumW{0};

    protected:
        //! \brief Neighborhood filter
        NeighborFilter m_nFilter;

        //! \brief Represent the current state of the fit (finalize function
        //! update the state)
        FIT_RESULT m_eCurrentState{UNDEFINED};

    public:
        /**************************************************************************/
        /* Initialization                                                         */
        /**************************************************************************/
        PONCA_FITTING_APIDOC_SETWFUNC
        PONCA_MULTIARCH inline void setNeighborFilter(const NeighborFilter& _nFilter) { m_nFilter = _nFilter; }

        PONCA_FITTING_APIDOC_INIT
        PONCA_MULTIARCH inline void init()
        {
            m_eCurrentState = UNDEFINED;
            startNewPass();
        }

        /*!
         * \brief Is the primitive well-fitted and ready to use (finalize has been called) ?
         * \warning The fit can still be unstable \see isStable
         */
        PONCA_MULTIARCH [[nodiscard]] inline bool isReady() const
        {
            return (m_eCurrentState == STABLE) || (m_eCurrentState == UNSTABLE);
        }

        /*! \brief Is the fitted primitive ready to use (finalize has been called and the result is stable) */
        PONCA_MULTIARCH [[nodiscard]] inline bool isStable() const { return m_eCurrentState == STABLE; }

        /*! \brief Is another pass required for fitting (finalize has been called and the result is #NEED_OTHER_PASS)
         * \see startNewPass */
        PONCA_MULTIARCH [[nodiscard]] inline bool needAnotherPass() const { return m_eCurrentState == NEED_OTHER_PASS; }

        /*! \brief Get number of points added in the neighborhood (with positive weight)  */
        PONCA_MULTIARCH [[nodiscard]] inline int getNumNeighbors() const { return m_nbNeighbors; }

        /*! \brief Get the sum of the weights */
        PONCA_MULTIARCH [[nodiscard]] inline Scalar getWeightSum() const { return m_sumW; }

        /*! \brief To be called when starting a new processing pass, ie. when `getCurrentState()==#NEED_ANOTHER_PASS` */
        PONCA_MULTIARCH inline void startNewPass()
        {
            m_nbNeighbors = 0;
            m_sumW        = Scalar(0);
        }

        /*! \brief Read access to the NeighborFilter \see setNeighborFilter */
        PONCA_MULTIARCH inline const NeighborFilter& getNeighborFilter() const { return m_nFilter; }

    public:
        /*! \return the current test of the fit */
        PONCA_MULTIARCH [[nodiscard]] inline FIT_RESULT getCurrentState() const { return m_eCurrentState; }

        PONCA_FITTING_APIDOC_ADDNEIGHBOR
        PONCA_MULTIARCH inline void addLocalNeighbor(Scalar w, const VectorType&, const DataPoint&)
        {
            m_sumW += w;
            ++(m_nbNeighbors);
        }

        PONCA_FITTING_APIDOC_FINALIZE
        PONCA_MULTIARCH inline FIT_RESULT finalize()
        {
            // handle specific configurations
            if (m_sumW == Scalar(0.) || m_nbNeighbors < 1)
            {
                return m_eCurrentState = UNDEFINED;
            }
            return m_eCurrentState = STABLE;
        }

    }; // class BasketUnit

    /**
        \brief Base class of any computation unit of a BasketDiff

        The differentiation can be done automatically in scale and/or space, by combining the enum values FitScaleDer
        and FitSpaceDer in the template parameter Type.

        The differentiated values are stored in static arrays. The size of the
        arrays is computed with respect to the derivation type (scale and/or space)
        and the number of the dimension of the ambiant space.
        By convention, the scale derivatives are stored at index 0 when Type
        contains at least FitScaleDer. The size of these arrays can be known using
        derDimension(), and the differentiation type by isScaleDer() and
        isSpaceDer().

        Thanks to the BasketDiff definition, we know that BasketDiffUnitBase has BasketUnitBase
        as base class (through the Basket). As a result, this class first asks to
        compute the Fit, and if it works properly, compute the weight derivatives.
     */
    template <class DataPoint, class _NFilter, int Type, typename T>
        requires ProvidesBasketUnitBase<T>
    class BasketDiffUnitBase : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_MATRIX_TYPE
    protected:
        static constexpr int NbDerivatives =
            ((Type & FitScaleDer) ? 1 : 0) + ((Type & FitSpaceDer) ? DataPoint::Dim : 0);
        static constexpr int DerStorageOrder = (Type & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor;

    public:
        /*! \brief Static array of scalars with a size adapted to the differentiation type */
        using VectorArray = Eigen::Matrix<Scalar, DataPoint::Dim, NbDerivatives, DerStorageOrder>;

        /*! \brief Static array of scalars with a size adapted to the differentiation type */
        using ScalarArray = Eigen::Matrix<Scalar, 1, NbDerivatives>;

    protected:
        // computation data
        ScalarArray m_dSumW; /*!< \brief Sum of weight derivatives */

    public:
        /************************************************************************/
        /* Initialization                                                       */
        /************************************************************************/
        /*! \see Concept::FittingProcedureConcept::init() */
        PONCA_MULTIARCH inline void init()
        {
            Base::init();
            m_dSumW.setZero();
        }

        /************************************************************************/
        /* Processing                                                           */
        /************************************************************************/
        /*! \see Concept::FittingProcedureConcept::addLocalNeighbor()
         * \todo update doc (takes derivatives) */
        PONCA_MULTIARCH inline void addLocalNeighbor(Scalar w, const VectorType& localQ, const DataPoint& attributes,
                                                     ScalarArray& dw)
        {
            Base::addLocalNeighbor(w, localQ, attributes);
            // call the Primitive Fit (without dw)
            int spaceId = (Type & FitScaleDer) ? 1 : 0;
            // compute weight
            if (Type & FitScaleDer)
                dw[0] = Base::getNeighborFilter().scaledw(attributes.pos(), attributes);

            if (Type & FitSpaceDer)
                dw.template segment<int(DataPoint::Dim)>(spaceId) =
                    -Base::getNeighborFilter().spacedw(attributes.pos(), attributes).transpose();

            m_dSumW += dw;
        }

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
        /*! \brief State specified at compilation time to compute derivatives in scale */
        PONCA_MULTIARCH [[nodiscard]] static inline constexpr bool isScaleDer() { return bool(Type & FitScaleDer); }
        /*! \brief State specified at compilation time to compute derivatives in space */
        PONCA_MULTIARCH [[nodiscard]] static inline constexpr bool isSpaceDer() { return bool(Type & FitSpaceDer); }
        /*! \brief Number of dimensions used for the differentiation */
        PONCA_MULTIARCH [[nodiscard]] static inline constexpr unsigned int derDimension() { return NbDerivatives; }
    };

} // namespace Ponca
