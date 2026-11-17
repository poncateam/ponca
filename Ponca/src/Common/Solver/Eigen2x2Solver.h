/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>

#include "../defines.h"

namespace Ponca
{
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

            PONCA_MULTIARCH [[nodiscard]] inline bool operator==(const Eigen2x2Solver<Scalar>& _other) const
            {
                return m_eivals.isApprox(_other.m_eivals) && m_eivecs.isApprox(_other.m_eivecs);
            }

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
    } // namespace internal
} // namespace Ponca

