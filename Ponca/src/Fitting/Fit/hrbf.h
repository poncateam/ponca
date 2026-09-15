/* This Source Code Form is subject to the terms of the Mozilla Public License,
 * v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at http://mozilla.org/MPL/2.0/.
 */
/// @file hrbf.h
/// This file is an adaptation of the HRBF-C++-Code by R. Vaillant, available at
/// https://rodolphe-vaillant.fr/entry/12/recipe-for-implicit-surface-reconstruction-with-hrbf
///
/// Original author: Gael Guennebaud - gael.guennebaud@inria.fr - http://www.labri.fr/perso/guenneba/
/// Rodolphe Vaillant - (Fixed the gradient evaluation) - http://www.irit.fr/~Rodolphe.Vaillant
///
/// Adaptation in Ponca formalism: Nicolas Mellado nicolas.mellado@irit.fr

#pragma once

#include "../defines.h"
#include "../concepts.h"
#include "../compute.h"
#include "../Filters/weightKernel.h"

#include <Eigen/Dense>

#define HRBF_REQUIREMENTS ProvidesNormal<P>

namespace Ponca
{

    /*!
     * \brief Corrected Normal Current Fit type.
     *
     * This fitting method generates triangles from a set a points cloud and use a statistical formula to compute :
     * - The principal curvatures values and directions
     * - The mean and gaussian curvatures
     *
     * \see ProvidesPrincipalCurvatures
     */
    template <class P, typename _WeightKernel = Pow3WeightKernel<typename P::Scalar>>
        requires HRBF_REQUIREMENTS
    class HRBF : ComputeObject<HRBF<P, _NeighborFilter>>
    {
    public:
        using DataPoint    = P;
        using MatrixType   = typename DataPoint::MatrixType;
        using MatrixDX     = typename Eigen::Matrix<Scalar, DataPoint::Dim, Eigen::Dynamic>;
        using Scalar       = typename DataPoint::Scalar;
        using VectorType   = typename DataPoint::VectorType;
        using DenseVector  = Eigen::VectorXd;
        using DenseMatrix  = Eigen::MatrixXd;
        using WeightKernel = _WeightKernel;

    protected:
        // Interpolation kernel
        WeightKernel m_wKernel;

        /// Each column represents p_i:  VectorX pi = _node_centers.col(i);
        MatrixDX m_node_centers;
        /// Vector of scalar values alpha
        DenseVector m_alphas;
        /// Each column represents beta_i:  VectorX bi = _betas.col(i);
        MatrixDX m_betas;

        FIT_RESULT m_eCurrentState{UNDEFINED};

    public:
        //! \brief Set the scalar field values to 0 and reset the isNormalized() status
        PONCA_MULTIARCH inline void init()
        {
            m_eCurrentState = UNDEFINED;
            m_node_centers.clear();
            m_alphas.clear();
            m_betas.clear();
        }

        /*!
         * \brief Convenience function for STL-like iterators
         * Add neighbors stored in a container using STL-like iterators, and call finalize at the end.
         * The fit is evaluated multiple time if needed (see #NEED_OTHER_PASS)
         */
        template <typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end)
        {
            constexpr int DIM       = DataPoint::Dim;
            int nb_points           = std::distance(begin, end);
            int nb_hrbf_constraints = (DIM + 1) * nb_points;
            int nb_constraints      = nb_hrbf_constraints;
            int nb_coeffs           = (DIM + 1) * nb_points;

            m_node_centers.resize(DIM, nb_points);
            m_betas.resize(DIM, nb_points);
            m_alphas.resize(nb_points);

            // Assemble the "design" and "value" matrix and vector
            DenseMatrix D(nb_constraints, nb_coeffs);
            DenseVector f(nb_constraints);
            DenseVector x(nb_coeffs);

            WeightKernel wk;

            // copy the node centers
            auto it = begin;
            for (int i = 0; i < nb_points; ++i, ++it)
                _node_centers.col(i) = (*it);

            it = begin;
            for (int i = 0; i < nb_points; ++i, ++it)
            {
                Vector p = (*it).pos();
                Vector n = (*it).normal();

                int io                          = (DIM + 1) * i;
                f(io)                           = 0;
                f.template segment<DIM>(io + 1) = n;

                for (int j = 0; j < nb_points; ++j)
                {
                    int jo          = (DIM + 1) * j;
                    VectorType diff = p - m_node_centers.col(j);
                    Scalar l        = diff.norm();
                    if (l == 0)
                    {
                        D.template block<DIM + 1, DIM + 1>(io, jo).setZero();
                    }
                    else
                    {
                        Scalar w                                   = wk.f(l);
                        Scalar dw_l                                = wk.df(l) / l;
                        Scalar ddw                                 = wk.ddf(l);
                        VectorType g                               = diff * dw_l;
                        D(io, jo)                                  = w;
                        D.row(io).template segment<DIM>(jo + 1)    = g.transpose();
                        D.col(jo).template segment<DIM>(io + 1)    = g;
                        D.template block<DIM, DIM>(io + 1, jo + 1) = (ddw - dw_l) / (l * l) * (diff * diff.transpose());
                        D.template block<DIM, DIM>(io + 1, jo + 1).diagonal().array() += dw_l;
                    }
                }
            }

            x = D.lu().solve(f);
            Eigen::Map<Eigen::Matrix<Scalar, DIM + 1, Eigen::Dynamic>> mx(x.data(), DIM + 1, nb_points);

            m_alphas = mx.row(0);
            m_betas  = mx.template bottomRows<DIM>();

            return m_eCurrentState = STABLE;
        }

        /*!
         * \brief Compute function that iterates over a subset of sampled points from an STL-Like container.
         * \tparam IndexRange An STL-like container storing the indices of the neighbors
         * \tparam PointContainer An STL-like container storing the points
         */
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(const IndexRange& ids, const PointContainer& points);

        PONCA_MULTIARCH [[nodiscard]] inline Scalar potential(const VectorType& x) const
        {
            Scalar ret   = 0;
            int nb_nodes = m_node_centers.cols();

            WeightKernel wk;

            for (int i = 0; i < nb_nodes; ++i)
            {
                VectorType diff = x - m_node_centers.col(i);
                Scalar l        = diff.norm();

                if (l > 0)
                {
                    ret += m_alphas(i) * wk.f(l);
                    ret += m_betas.col(i).dot(diff) * wk.df(l) / l;
                }
            }
            return ret;
        }

        VectorType primitiveGradient(const Vector& x) const
        {
            VectorType grad = VectorType::Zero();
            int nb_nodes    = m_node_centers.cols();

            WeightKernel wk;

            for (int i = 0; i < nb_nodes; i++)
            {
                VectorType node = m_node_centers.col(i);
                VectorType beta = m_betas.col(i);
                Scalar alpha    = m_alphas(i);
                VectorType diff = x - node;

                VectorType diffNormalized = diff;
                Scalar l                  = diff.norm();

                if (l > 0.00001f)
                {
                    diffNormalized.normalize();
                    Scalar dphi  = wk.df(l);
                    Scalar ddphi = wk.ddf(l);

                    Scalar alpha_dphi = alpha * dphi;

                    Scalar bDotd_l   = beta.dot(diff) / l;
                    Scalar squared_l = diff.squaredNorm();

                    grad += alpha_dphi * diffNormalized;
                    grad += bDotd_l * (ddphi * diffNormalized - diff * dphi / squared_l) + beta * dphi / l;
                }
            }
            return grad;
        }

        //
        // //! \brief Comparison operator
        // PONCA_MULTIARCH [[nodiscard]] bool operator==(const CNC& other) const
        // {
        //     // We use the matrix to compare the fitting results
        //     return (m_eCurrentState == other.m_eCurrentState) && ;
        // }

        // //! \brief Comparison operator, convenience function
        // PONCA_MULTIARCH [[nodiscard]] bool operator!=(const CNC& other) const
        // {
        //     // We use the matrix to compare the fitting results
        //     return !(this == &other);
        // }

        // //! \brief Approximate operator
        // PONCA_MULTIARCH [[nodiscard]] bool isApprox(
        //     const CNC& other, const Scalar& epsilon = Eigen::NumTraits<Scalar>::dummy_precision()) const
        // {
        //     PONCA_MULTIARCH_STD_MATH(abs);
        //
        //     return (m_eCurrentState == other.m_eCurrentState) && (std::abs(kMean() - other.kMean()) < epsilon) &&
        //            (std::abs(GaussianCurvature() - other.GaussianCurvature()) < epsilon) &&
        //            (std::abs(kmin() - other.kmin()) < epsilon) && (std::abs(kmax() - other.kmax()) < epsilon);
        // }

        //! \brief Is the fitted primitive ready to use (finalize has been called and the result is stable)
        PONCA_MULTIARCH [[nodiscard]] inline bool isStable() const { return m_eCurrentState == STABLE; }

        // TODO Add methods to respect concept ProvidesImplicitPrimitive
    }; // class CNC

} // namespace Ponca

#include "cnc.hpp"
