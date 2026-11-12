/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../defines.h"
#include "../concepts.h"

#define FUNDAMENTAL_FORM_WEINGARTEN_ESTIMATOR_REQUIREMENTS \
    ProvidesFirstFondamentalFormComponents<T>&& ProvidesSecondFondamentalFormComponents<T>&& Is3D<DataPoint>
#define NORMAL_DERIVATIVE_WEINGARTEN_ESTIMATOR_REQUIREMENTS Is3D<DataPoint>&& ProvidesSpaceDerivative<T>
#define WIENGARTEN_CURVATURE_ESTIMATOR_REQUIREMENTS \
    ProvidesTangentPlaneBasis<T>&& ProvidesWeingartenMap<T>&& Is3D<DataPoint>

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
        \verbatim ProvidesFirstFondamentalFormComponents, ProvidesSecondFondamentalFormComponents \endverbatim
        */
    template <class DataPoint, class _NFilter, typename T>
        requires FUNDAMENTAL_FORM_WEINGARTEN_ESTIMATOR_REQUIREMENTS
    class FundamentalFormWeingartenEstimator : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;

    public:
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS(FundamentalFormWeingartenEstimator, fundamentalFormWeingartenEstimator)
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS(FundamentalFormWeingartenEstimator, firstFondamentalFormComponent)
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS(FundamentalFormWeingartenEstimator, secondFondamentalFormComponent)

        /// \brief Assembles and returns the first fundamental form from the base class
        ///
        /// \return first fundamental form
        /// \see firstFundamentalForm(Matrix2& first) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 firstFundamentalForm() const;

        /// \copybrief firstFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template <typename Matrix2Derived>
        PONCA_MULTIARCH inline void firstFundamentalForm(Matrix2Derived& first) const;

        /// \brief Assembles and returns the second fundamental form from the base class
        ///
        /// \return second fundamental form
        /// \see secondFundamentalForm(Matrix2& second) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 secondFundamentalForm() const;

        /// \copybrief secondFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template <typename Matrix2Derived>
        PONCA_MULTIARCH inline void secondFundamentalForm(Matrix2Derived& second) const;

        //! \brief Returns the Weingarten Map
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template <typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived& w) const;

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

        \note Computations are marked as `UNSTABLE` if the computed basis does not properly align with the gradient of
       the fitted primitive.

        This primitive provides:
        \verbatim PROVIDES_WEINGARTEN_MAP, ProvidesTangentPlaneBasis\endverbatim

        This primitive requires:
        \verbatim ProvidesNormalDerivative \endverbatim
        */
    template <class DataPoint, class _NFilter, int DiffType, typename T>
        requires NORMAL_DERIVATIVE_WEINGARTEN_ESTIMATOR_REQUIREMENTS
    class NormalDerivativeWeingartenEstimator : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_MATRIX_TYPE
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;

    private:
        MatrixType m_tangentBasis{MatrixType::Zero()};

    public:
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS_DER(NormalDerivativeWeingartenEstimator,
                                                 normalDerivativeWeingartenEstimator)
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS_DER(NormalDerivativeWeingartenEstimator, tangentPlaneBasis)
        PONCA_FITTING_DECLARE_FINALIZE

        //! \brief Returns the Weingarten Map
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template <typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived& w) const;

        /// \copydoc CovariancePlaneFitImpl::worldToTangentPlane
        PONCA_MULTIARCH inline VectorType worldToTangentPlane(const VectorType& _q,
                                                              bool _isPositionVector = true) const;

        /// \copydoc CovariancePlaneFitImpl::tangentPlaneToWorld
        PONCA_MULTIARCH inline VectorType tangentPlaneToWorld(const VectorType& _q,
                                                              bool _isPositionVector = true) const;
    };

    namespace internal
    {
        /*!
            \brief Closed-form eigen decomposition of a 2x2 matrix with real eigenvalues, not necessarily symmetric

            Drop-in replacement for `Eigen::SelfAdjointEigenSolver<Matrix2>::computeDirect` but for cases where the 
            matrix is not symmetric. Unlike other eigen slovers, this code works within cuda kernels.  
        */
        template <typename Scalar>
        class Eigen2x2Solver
        {
        public:
            using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
            using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

            PONCA_MULTIARCH inline void computeDirect(const Matrix2& _w)
            {
                PONCA_MULTIARCH_STD_MATH(sqrt);
                PONCA_MULTIARCH_STD_MATH(max);

                // Eigenvalues: roots of k^2 - tr(w) k + det(w) = 0
                const Scalar tr    = Scalar(0.5) * (_w(0, 0) + _w(1, 1)); // 0.5 here reduces mult count
                const Scalar det   = _w(0, 0) * _w(1, 1) - _w(0, 1) * _w(1, 0);
                const Scalar delta = sqrt(max(tr * tr - det, Scalar(0)));
                m_eivals << tr - delta, tr + delta;

                const Scalar tol = Eigen::NumTraits<Scalar>::dummy_precision() * (_w.cwiseAbs().maxCoeff() + Scalar(1));

                // Eigenvectors, pick an orthonormal pair when w = kI as every vector is an eigenvector...                
                // (v1 first, so that the fallback for v0 can be made orthogonal to it)
                Vector2 v0, v1;
                if (!eigenvector(_w, m_eivals(1), tol, v1))
                    v1 = Vector2::UnitY();
                // A symmetric matrix has orthogonal eigenvectors: take v0 orthogonal to v1 instead of computing it, as
                // computing it independently is ill-conditioned near umbilics (eigengap ~ 0)
                if (_w(0, 1) == _w(1, 0) || !eigenvector(_w, m_eivals(0), tol, v0))
                    v0 = Vector2(v1.y(), -v1.x());

                m_eivecs.col(0) = v0;
                m_eivecs.col(1) = v1;
            }

            //! \brief Eigenvalues, in increasing order
            PONCA_MULTIARCH [[nodiscard]] inline const Vector2& eigenvalues() const { return m_eivals; }

            //! \brief Normalized eigenvectors (columns), in the order of eigenvalues()
            PONCA_MULTIARCH [[nodiscard]] inline const Matrix2& eigenvectors() const { return m_eivecs; }

        private:
            Vector2 m_eivals{Vector2::Zero()};
            Matrix2 m_eivecs{Matrix2::Identity()};

            //! \brief Normalized eigenvector of _w for its eigenvalue _k
            PONCA_MULTIARCH static inline bool eigenvector(const Matrix2& _w, Scalar _k, Scalar _tol, Vector2& _v)
            {
                const Vector2 v0(_w(0, 1), _k - _w(0, 0)); // 1st candidate, 0 if w = kI
                const Vector2 v1(_k - _w(1, 1), _w(1, 0)); // 2cd candidate, 0 if w = kI

                // Use the one with the highest norm avoid numerical instabilities
                _v = v0.squaredNorm() > v1.squaredNorm() ? v0 : v1;

                const Scalar n = _v.norm();
                if (n <= _tol)
                    return false;

                _v /= n;
                return true;
            }
        };

        /*!
            \brief Compute principal curvatures from a base class providing fundamental forms

            Since the weingarten map is self-adjoint, a basis exists consisting of the eigenvectors of W.
            The eigenvectors, combined with their eigenvalues, define the principal curvature directions and
            values.
            The Gaussian and mean curvatures can be retrieved as \f$\text{K} = \text{det}(\text{W})\f$
            and \f$\text{H} = \frac{1}{2}\text{trace}(\text{W})\f$, respectively.


            \warning To map from the tangent plane to ambiant space, we assume that the patch is represented as
            (h(u,v), u, v) because CovariancePlaneFitImpl::worldToTangentPlane() wraps up coordinates in this order
           (height, u and v).

            This primitive requires:
            \verbatim ProvidesTangentPlaneBasis, PROVIDES_WEINGARTEN_MAP, ProvidesPrincipalCurvatures \endverbatim
            */
        template <class DataPoint, class _NFilter, typename T>
            requires WIENGARTEN_CURVATURE_ESTIMATOR_REQUIREMENTS
        class WeingartenCurvatureEstimatorBase : public T
        {
            PONCA_FITTING_DECLARE_DEFAULT_TYPES
            using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;

        private:
            Eigen2x2Solver<Scalar> m_solver; // W is not necessarily symmetric, see Eigen2x2Solver

        public:
            PONCA_FITTING_DECLARE_FINALIZE

            //! \brief Returns an estimate of the minimal principal curvature value
            PONCA_MULTIARCH [[nodiscard]] inline Scalar kmin() const { return m_solver.eigenvalues().x(); }

            //! \brief Returns an estimate of the maximal principal curvature value
            PONCA_MULTIARCH [[nodiscard]] inline Scalar kmax() const { return m_solver.eigenvalues().y(); }

            //! \brief Returns an estimate of the minimal principal curvature direction
            PONCA_MULTIARCH [[nodiscard]] inline VectorType kminDirection() const
            {
                VectorType vmin;
                vmin(0)                       = Scalar(0); // set height
                vmin.template bottomRows<2>() = m_solver.eigenvectors().col(0);
                vmin                          = Base::tangentPlaneToWorld(vmin, false);
                return vmin;
            }

            //! \brief Returns an estimate of the maximal principal curvature direction
            PONCA_MULTIARCH [[nodiscard]] inline VectorType kmaxDirection() const
            {
                VectorType vmax;
                vmax(0)                       = Scalar(0); // set height
                vmax.template bottomRows<2>() = m_solver.eigenvectors().col(1);
                vmax                          = Base::tangentPlaneToWorld(vmax, false);
                return vmax;
            }

            //! \brief Returns an estimate of the mean curvature
            PONCA_MULTIARCH [[nodiscard]] inline Scalar kMean() const { return (kmin() + kmax()) / Scalar(2); }

            //! \brief Returns an estimate of the Gaussian curvature
            PONCA_MULTIARCH [[nodiscard]] inline Scalar GaussianCurvature() const { return kmin() * kmax(); }
        };
    } // namespace internal

    /// \copydoc internal::WeingartenCurvatureEstimatorBase
    template <class DataPoint, class _NFilter, typename T>
    struct WeingartenCurvatureEstimator : public internal::WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS(WeingartenCurvatureEstimator, meanCurvature)
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS(WeingartenCurvatureEstimator, curvatureTensor)
    };

    /// \copydoc internal::WeingartenCurvatureEstimatorBase
    template <class DataPoint, class _NFilter, int DiffType, typename T>
    struct WeingartenCurvatureEstimatorDer : public internal::WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS_DER(WeingartenCurvatureEstimatorDer, meanCurvature)
        PONCA_EXPLICIT_BASKET_CAST_OPERATORS_DER(WeingartenCurvatureEstimatorDer, curvatureTensor)
    };

} // namespace Ponca

#include "weingarten.hpp"

