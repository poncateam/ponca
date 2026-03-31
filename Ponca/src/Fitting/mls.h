/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "defines.h"
#include "enums.h"

#include "./compute.h"

#include <algorithm>

namespace Ponca
{
    /*!
     * \brief Computes the fit using the MLS iteration process \cite oztireli:2009:feature
     *
     * The position of the projected point is outputted within getNeighborFilter().evalPos()
     *
     * \tparam Scalar scalar type
     */
    template <typename Scalar>
    struct MLS
    {
    public:
        MLS()
        {
            setPrecision(Eigen::NumTraits<Scalar>::dummy_precision());
            setNIter(5);
        }

        MLS(unsigned int nIter) : MLS() { setNIter(nIter); }

        MLS(unsigned int nIter, Scalar eps) : MLS(nIter) { setPrecision(Eigen::NumTraits<Scalar>::dummy_precision()); }

        void setPrecision(Scalar newEps)
        {
            // Enforce positiveness
            m_eps = std::abs(newEps);
        }

        void setNIter(unsigned int m_nIter)
        {
            // At least 1, otherwise no computation would be performed
            m_nIter = std::max(m_nIter, 1u);
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam ItB Begin iterator type
         * \tparam ItE End iterator type
         *
         * \param fit The fitting object
         * \param beg The begining of point range
         * \param end The end of point range
         *
         * \return The result of the fit
         */
        template <typename Fit, typename ItB, typename ItE>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& fit, const ItB& itb, const ItE& ite)
        {
            return computeMLSImpl(fit, [&]() { return compute(fit, itb, ite); });
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam Container Container of points
         *
         * \param fit The fitting object
         * \param container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& fit, const PointContainer& container)
        {
            return compute(fit, std::begin(container), std::end(container));
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam IdxRange Range index type
         * \tparam PointContainer Point container (must provide random access)
         *
         * \param fit The fitting object
         * \param range The container of indices
         * \param container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename IdxRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(Fit& fit, const IdxRange& range,
                                                         const PointContainer& container)
        {
            return computeMLSImpl(fit, [&]() { return computeWithIds(fit, range, container); });
        }

    protected:
        Scalar m_eps         = Eigen::NumTraits<Scalar>::dummy_precision();
        unsigned int m_nIter = 5;

        /*!
         * \brief Computes the fit using the MLS iteration process.
         *
         * The position of the projected point is outputted within getNeighborFilter().evalPos()
         *
         * \tparam Fit The fitting type
         * \tparam Func The compute procedure
         *
         * \param fit The fitting object
         * \param compute The procedure to compute estimator
         *
         * \return The result of the fit
         */
        template <typename Fit, typename Func>
        PONCA_MULTIARCH inline FIT_RESULT computeMLSImpl(Fit& fit, Func&& compute)
        {
            FIT_RESULT res = UNDEFINED;
            auto filter    = fit.getNeighborFilter();
            auto lastPos   = filter.evalPos();

            for (unsigned int mm = 0; mm < m_nIter; ++mm)
            {
                filter.changeNeighborhoodFrame(lastPos);
                fit.setNeighborFilter(filter);
                res = compute();

                if (fit.isStable())
                {
                    auto newPos = fit.project(lastPos);
                    if (newPos.isApprox(lastPos, m_eps))
                        return res;
                    lastPos = newPos;
                }
                else
                {
                    return res;
                }
            }

            return res;
        }
    };
} // namespace Ponca
