/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Core>

namespace Ponca
{
    template < class DataPoint, class _NFilter, typename T>
/**
 *
 * \brief Base class for any 3d curvature estimator: holds \f$k_{\min}\f$, \f$k_{\max}\f$ and associated vectors,
 * such that \f$ k_{\min} <= k_{\max} \f$
 */
    class CurvatureEstimatorBase : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

    protected:
        enum
        {
            PROVIDES_PRINCIPAL_CURVATURES
        };

    private:
        /// \brief Minimal principal curvature
        Scalar m_kmin {0},
        /// \brief Maximal principal curvature
        m_kmax {0};
        /// \brief Direction associated to the minimal principal curvature
        VectorType m_vmin {VectorType::Zero()},
        /// \brief Direction associated to the maximal principal curvature
        m_vmax {VectorType::Zero()};

        /// \brief Internal state indicating if the curvature are set to default (false), or have been computed (true)
        bool m_isValid {false};

        static_assert ( DataPoint::Dim == 3, "CurvatureEstimatorBase is only valid in 3D");

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(CurvatureEstimatorBase, curvatureEstimatorBase)
        PONCA_FITTING_DECLARE_INIT

        /// \brief Returns true if contains valid curvature values (and not default ones)
        PONCA_MULTIARCH [[nodiscard]] inline bool isValid() const { return m_isValid; }

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
        //! \brief Returns an estimate of the minimal principal curvature value
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kmin() const { return m_kmin; }

        //! \brief Returns an estimate of the maximal principal curvature value
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kmax() const { return m_kmax; }

        //! \brief Returns an estimate of the minimal principal curvature direction
        PONCA_MULTIARCH [[nodiscard]] inline VectorType kminDirection() const { return m_vmin; }

        //! \brief Returns an estimate of the maximal principal curvature direction
        PONCA_MULTIARCH [[nodiscard]] inline VectorType kmaxDirection() const { return m_vmax; }

        //! \brief Returns an estimate of the mean curvature
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kMean() const { return (m_kmin + m_kmax)/Scalar(2);}

        //! \brief Returns an estimate of the Gaussian curvature
        PONCA_MULTIARCH [[nodiscard]] inline Scalar GaussianCurvature() const { return m_kmin * m_kmax;}

    protected:
        /// \brief Set curvature values. To be called in finalize() by child classes
        ///
        /// If the given parameters are such that \f$ k_{\min} > k_{\max} \f$, then this
        /// method swap the two curvature values and directions and store them such that
        /// \f$ k_{\min} <= k_{\max} \f$.
        ///
        PONCA_MULTIARCH inline void setCurvatureValues(Scalar kmin, Scalar kmax, const VectorType& vmin, const VectorType& vmax);
    };

    template < class DataPoint, class WFunctor, int /*DiffType*/, typename T>
    using CurvatureEstimatorBaseDiff = CurvatureEstimatorBase<DataPoint, WFunctor, T>;



/*!
    \brief Compute principal curvatures from a base class providing fundamental forms

    \see MongePatchPrimitive

    This primitive provides:
    \verbatim PROVIDES_PRINCIPAL_CURVATURES \endverbatim

    This primitive requires:
    \verbatim PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class FundamentalFormCurvatureEstimator : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "FundamentalFormCurvatureEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS &&
                    Base::PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS,
            PROVIDES_PRINCIPAL_CURVATURES  /*!< \brief Provides curvature API */
        };

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(FundamentalFormCurvatureEstimator, fundamentalFormCurvatureEstimator)

        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 firstFundamentalForm() const {
            Matrix2 first;
            Base::firstFundamentalFormComponents(*(first.data()), *(first.data()+1), *(first.data()+3));
            *(first.data()+2) = *(first.data()+1);
            return first;
        }

        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 secondFundamentalForm() const {
            Matrix2 second;
            Base::secondFundamentalFormComponents(*(second.data()), *(second.data()+1), *(second.data()+3));
            *(second.data()+2) = *(second.data()+1);
            return second;
        }

        //! \brief Returns the WeingartenMap, defined as
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const {
            return firstFundamentalForm().inverse() *  secondFundamentalForm();
        }

//        //! \brief Returns an estimate of the minimal principal curvature value
//        PONCA_MULTIARCH inline Scalar kmin() const { return m_kmin; }
//
//        //! \brief Returns an estimate of the maximal principal curvature value
//        PONCA_MULTIARCH inline Scalar kmax() const { return m_kmax; }
//
//        //! \brief Returns an estimate of the minimal principal curvature direction
//        PONCA_MULTIARCH inline VectorType kminDirection() const { return m_vmin; }
//
//        //! \brief Returns an estimate of the maximal principal curvature direction
//        PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_vmax; }

        //! \brief Returns an estimate of the mean curvature
        PONCA_MULTIARCH inline Scalar kMean() const {
            Scalar E, F, G, L, M, N;
            Base::firstFundamentalFormComponents(E, F, G);
            Base::secondFundamentalFormComponents(L, M, N);
            return (G*L-Scalar(2)*F*M+E*N)/(Scalar(2)*(E*G-F*F));
//            return Scalar(0.5)*weingartenMap().trace();
        }

        //! \brief Returns an estimate of the Gaussian curvature
        PONCA_MULTIARCH inline Scalar GaussianCurvature() const {
            Scalar E, F, G, L, M, N;
            Base::firstFundamentalFormComponents(E, F, G);
            Base::secondFundamentalFormComponents(L, M, N);
            return (L*N-M*M)/(E*G-F*F);
//            return weingartenMap().determinant();
        }
    };

    template < class DataPoint, class _NFilter, int /*DiffType*/, typename T>
    using FundamentalFormCurvatureEstimatorDiff = FundamentalFormCurvatureEstimator<DataPoint, _NFilter, T>;
} //namespace Ponca

#include "curvature.hpp"
