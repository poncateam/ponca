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
#include "../Filters/weightFilter.h"

#include <Eigen/Dense>

#define HRBF_REQUIREMENTS ProvidesNormal<P>

namespace Ponca
{

    /*!
     * \brief Hermite Radial Basis Function (HRBF) implicit surface interpolation.
     *
     * Builds a scalar field \f$ f \f$ that interpolates both the positions and the normals of a set of oriented
     * points.
     *
     * \tparam P Point type, must provide normals (see ProvidesNormal)
     * \tparam _WeightKernel Interpolation kernel. Defaults to Pow3WeightKernel.
     */
    template <class P, typename _WeightKernel = Pow3WeightKernel<typename P::Scalar>>
        requires HRBF_REQUIREMENTS
    class HRBF : ComputeObject<HRBF<P, _WeightKernel>>
    {
    public:
        using DataPoint    = P;
        using Scalar       = typename DataPoint::Scalar;
        using MatrixType   = typename DataPoint::MatrixType;
        using VectorType   = typename DataPoint::VectorType;
        using MatrixDX     = typename Eigen::Matrix<Scalar, DataPoint::Dim, Eigen::Dynamic>;
        using DenseVector  = Eigen::Vector<Scalar, Eigen::Dynamic>;
        using DenseMatrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using WeightKernel = _WeightKernel;

        static constexpr Scalar eps = Scalar(0.00001f);

    protected:
        /**
         * Very small class that emulates what hrbf needs from points. This class is used as
         * part as computeWithIds which requires acceptance of single pass iterator to work with 
         * kdtree, yet needs to compute data for each pairs of points. 
         */
        struct PointAdapter
        {
            VectorType _pos;
            VectorType _normal;

            PointAdapter(const VectorType& _p, const VectorType& _n) : _pos(_p), _normal(_n) {}

            const VectorType& pos() const { return _pos; }
            const VectorType& normal() const { return _normal; }
        };

        /// Each column represents p_i:  VectorX pi = _node_centers.col(i);
        MatrixDX m_node_centers;
        /// Vector of scalar values alpha
        DenseVector m_alphas;
        /// Each column represents beta_i:  VectorX bi = _betas.col(i);
        MatrixDX m_betas;

        FIT_RESULT m_eCurrentState{UNDEFINED};

    public:
        //! \brief Set the scalar field values to 0 and reset the isNormalized() status
        inline void init()
        {
            m_eCurrentState = UNDEFINED;

            // Clear matrices to make sure they can't be reused
            m_node_centers.resize(DataPoint::Dim, 0);
            m_alphas.resize(0);
            m_betas.resize(DataPoint::Dim, 0);
        }

        /*!
         * \brief Convenience function for STL-like iterators
         * Add neighbors stored in a container using STL-like iterators, and call finalize at the end.
         * The fit is evaluated multiple time if needed (see #NEED_OTHER_PASS)
         */
        template <typename IteratorBegin, typename IteratorEnd>
        inline FIT_RESULT compute(const IteratorBegin& begin, const IteratorEnd& end)
        {
            constexpr int DIM = DataPoint::Dim;

            const int nb_points           = std::distance(begin, end);
            const int nb_hrbf_constraints = (DIM + 1) * nb_points;
            const int nb_constraints      = nb_hrbf_constraints;
            const int nb_coeffs           = (DIM + 1) * nb_points;
            const WeightKernel wk;

            m_node_centers.resize(DIM, nb_points);
            m_betas.resize(DIM, nb_points);
            m_alphas.resize(nb_points);

            // Assemble the "design" and "value" matrix and vector
            DenseMatrix D(nb_constraints, nb_coeffs);
            DenseVector f(nb_constraints);

            // copy the node centers
            auto it = begin;
            for (int i = 0; i < nb_points; ++i, ++it)
                m_node_centers.col(i) = (*it).pos();

            it = begin;
            for (int i = 0; i < nb_points; ++i, ++it)
            {
                const auto& p = (*it).pos();
                const auto& n = (*it).normal();

                const int io                    = (DIM + 1) * i;
                f(io)                           = 0;
                f.template segment<DIM>(io + 1) = n;

                for (int j = 0; j < nb_points; ++j)
                {
                    const int jo          = (DIM + 1) * j;
                    const VectorType diff = p - m_node_centers.col(j);
                    const Scalar l        = diff.norm();

                    if (l < eps)
                    {
                        D.template block<DIM + 1, DIM + 1>(io, jo).setZero();
                    }
                    else
                    {
                        const Scalar w     = wk.f(l);
                        const Scalar dw_l  = wk.df(l) / l;
                        const Scalar ddw   = wk.ddf(l);
                        const VectorType g = diff * dw_l;

                        D(io, jo)                                  = w;
                        D.row(io).template segment<DIM>(jo + 1)    = g.transpose();
                        D.col(jo).template segment<DIM>(io + 1)    = g;
                        D.template block<DIM, DIM>(io + 1, jo + 1) = (ddw - dw_l) / (l * l) * (diff * diff.transpose());
                        D.template block<DIM, DIM>(io + 1, jo + 1).diagonal().array() += dw_l;
                    }
                }
            }

            const DenseVector x = D.lu().solve(f);
            Eigen::Map<const Eigen::Matrix<Scalar, DIM + 1, Eigen::Dynamic>> mx(x.data(), DIM + 1, nb_points);

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
        inline FIT_RESULT computeWithIds(const IndexRange& ids, const PointContainer& points)
        {
            std::vector<PointAdapter> adapter;

            auto begin = points.begin();
            for (const auto& idx : ids)
            {
                auto current = std::next(begin, idx);
                adapter.emplace_back(current->pos(), current->normal());
            }

            return compute(adapter.begin(), adapter.end());
        }

        [[nodiscard]] inline Scalar potential(const VectorType& x) const
        {
            const int nb_nodes = m_node_centers.cols();
            const WeightKernel wk;

            Scalar ret = 0;
            for (int i = 0; i < nb_nodes; ++i)
            {
                const VectorType diff = x - m_node_centers.col(i);
                const Scalar l        = diff.norm();

                if (l > eps)
                {
                    ret += m_alphas(i) * wk.f(l);
                    ret += m_betas.col(i).dot(diff) * wk.df(l) / l;
                }
            }
            return ret;
        }

        [[nodiscard]] inline VectorType primitiveGradient(const VectorType& x) const
        {
            const int nb_nodes = m_node_centers.cols();
            const WeightKernel wk;

            VectorType grad = VectorType::Zero();
            for (int i = 0; i < nb_nodes; i++)
            {
                const VectorType node = m_node_centers.col(i);
                const VectorType beta = m_betas.col(i);
                const VectorType diff = x - node;
                const Scalar alpha    = m_alphas(i);
                const Scalar l        = diff.norm();

                if (l > eps)
                {
                    const VectorType diffNormalized = diff / l;

                    const Scalar dphi  = wk.df(l);
                    const Scalar ddphi = wk.ddf(l);

                    const Scalar alpha_dphi = alpha * dphi;

                    const Scalar bDotd_l   = beta.dot(diff) / l;
                    const Scalar squared_l = diff.squaredNorm();

                    grad += alpha_dphi * diffNormalized;
                    grad += bDotd_l * (ddphi * diffNormalized - diff * dphi / squared_l) + beta * dphi / l;
                }
            }
            return grad;
        }

        //! \brief Is the fitted primitive ready to use (finalize has been called and the result is stable)
        PONCA_MULTIARCH [[nodiscard]] inline bool isStable() const { return m_eCurrentState == STABLE; }

        /// writing implementations of function that are required to use HRBF fitting alongside Basket)
        using NeighborFilter = DistWeightFilter<DataPoint, ConstantWeightKernel<Scalar>>;

        PONCA_FITTING_APIDOC_SETWFUNC
        PONCA_MULTIARCH inline bool setNeighborFilter(const NeighborFilter&) { return false; }

        PONCA_MULTIARCH inline const VectorType& project(const VectorType& x) const { return x; }

        PONCA_FITTING_IS_SIGNED(true)
    }; // class HRBF

} // namespace Ponca

