/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

namespace Ponca
{
/*!
    \brief Compute a Weingarten map from fundamental forms.

    The Weingarten map, also called Shape Operator, is defined in its matrix form as
    \f[ \text{W} = \text{F}_\text{I}^{-1}\text{F}_\text{II}, \f]
    with \f$\text{F}_\text{I}, \text{F}_\text{II}\f$ the first and second fundamental forms,
    respectively.

    \see WeingartenCurvatureEstimator and WeingartenCurvatureEstimatorDer to compute curvatures from the
    weingarten map

    This class also provides mean and Gaussian curvature directly from the fundamental forms.

    This primitive provides:
    \verbatim PROVIDES_WEINGARTEN_MAP \endverbatim

    This primitive requires:
    \verbatim PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class FundamentalFormWeingartenEstimator : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "FundamentalFormWeingartenEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS &&
                    Base::PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS,
            PROVIDES_WEINGARTEN_MAP
        };
    public:
        PONCA_EXPLICIT_CAST_OPERATORS(FundamentalFormWeingartenEstimator, fundamentalFormWeingartenEstimator)

        /// \brief Assembles and returns the first fundamental form from the base class
        ///
        /// \return first fundamental form
        /// \see firstFundamentalForm(Matrix2& first) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 firstFundamentalForm() const;

        /// \copybrief firstFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void firstFundamentalForm(Matrix2Derived& first) const;

        /// \brief Assembles and returns the second fundamental form from the base class
        ///
        /// \return second fundamental form
        /// \see secondFundamentalForm(Matrix2& second) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 secondFundamentalForm() const;

        /// \copybrief secondFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void secondFundamentalForm(Matrix2Derived &second) const;

        //! \brief Returns the Weingarten Map
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived &w) const;

        /// \brief Returns an estimate of the mean curvature directly from the fundamental forms
        PONCA_MULTIARCH inline Scalar kMean() const;

        //! \brief Returns an estimate of the Gaussian curvature directly from the fundamental forms
        PONCA_MULTIARCH inline Scalar GaussianCurvature() const;
    };


/*!
    \brief Compute a Weingarten map from the spatial derivatives of the normal field \f$ N \f$.

    \see WeingartenCurvatureEstimator and WeingartenCurvatureEstimatorDer to compute curvatures from the
        weingarten map

    A tangent plane aligned with the map is also computed during the process.

    \note Computations are marked as `UNSTABLE` if the computed basis does not properly align with the gradient of the
          fitted primitive.

    This primitive provides:
    \verbatim PROVIDES_WEINGARTEN_MAP, PROVIDES_TANGENT_PLANE_BASIS\endverbatim

    This primitive requires:
    \verbatim PROVIDES_NORMAL_DERIVATIVE \endverbatim
    */
    template < class DataPoint, class _NFilter, int DiffType, typename T >
    class NormalDerivativeWeingartenEstimator : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_MATRIX_TYPE
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "NormalDerivativeWeingartenEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_NORMAL_DERIVATIVE,
            PROVIDES_WEINGARTEN_MAP,
            PROVIDES_TANGENT_PLANE_BASIS
        };

    private:
        MatrixType m_tangentBasis {MatrixType::Zero()};

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(NormalDerivativeWeingartenEstimator, normalDerivativeWeingartenEstimator)
        PONCA_FITTING_DECLARE_FINALIZE

        //! \brief Returns the Weingarten Map
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived &w) const;


        /// \copydoc CovariancePlaneFitImpl::worldToTangentPlane
        PONCA_MULTIARCH inline VectorType worldToTangentPlane(const VectorType &_q,
                                                              bool _isPositionVector = true) const;

        /// \copydoc CovariancePlaneFitImpl::tangentPlaneToWorld
        PONCA_MULTIARCH inline VectorType tangentPlaneToWorld(const VectorType &_q,
                                                              bool _isPositionVector = true) const;
    };


namespace internal
{
    /*!
        \brief Compute principal curvatures from a base class providing fundamental forms

        Since the weingarten map is self-adjoint, a basis exists consisting of the eigenvectors of W.
        The eigenvectors, combined with their eigenvalues, define the principal curvature directions and
        values.
        The Gaussian and mean curvatures can be retrieved as \f$\text{K} = \text{det}(\text{W})\f$
        and \f$\text{H} = \frac{1}{2}\text{trace}(\text{W})\f$, respectively.


        \warning To map from the tangent plane to ambiant space, we assume that the patch is represented as
        (h(u,v), u, v) because CovariancePlaneFitImpl::worldToTangentPlane() wraps up coordinates in this order (height, u and v).

        This primitive requires:
        \verbatim PROVIDES_TANGENT_PLANE_BASIS, PROVIDES_WEINGARTEN_MAP, PROVIDES_PRINCIPAL_CURVATURES \endverbatim
        */
    template < class DataPoint, class _NFilter, typename T >
    class WeingartenCurvatureEstimatorBase : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
            using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "WeingartenCurvatureEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_TANGENT_PLANE_BASIS &&  // required for tangentPlaneToWorld
                    Base::PROVIDES_PRINCIPAL_CURVATURES && // required curvature storage
                    Base::PROVIDES_WEINGARTEN_MAP
        };

    public:
        PONCA_FITTING_DECLARE_FINALIZE
    };
}

    /// \copydoc internal::WeingartenCurvatureEstimatorBase
    template < class DataPoint, class _NFilter, typename T >
    struct WeingartenCurvatureEstimator : public internal::WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>{
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_EXPLICIT_CAST_OPERATORS(WeingartenCurvatureEstimator, weingartenCurvatureEstimator)
    };

    /// \copydoc internal::WeingartenCurvatureEstimatorBase
    template < class DataPoint, class _NFilter, int DiffType, typename T >
    struct WeingartenCurvatureEstimatorDer : public internal::WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>{
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        PONCA_EXPLICIT_CAST_OPERATORS_DER(WeingartenCurvatureEstimatorDer, weingartenCurvatureEstimator)

    };

} //namespace Ponca

#include "weingarten.hpp"

