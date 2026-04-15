/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "defines.h"
#include "enums.h"

#include "./compute.h"
#include "./project.h"

#include <algorithm>

namespace Ponca
{
    /*!
     * \brief Computes the fit using the Moving Least Squares approach.
     * The projection operator can be customized, see \cite Alexa:2004:projection.
     *
     * The position of the projected point is outputted within getNeighborFilter().evalPos()
     *
     * \tparam Scalar scalar type
     */
    template <typename Scalar>
    struct MLSEvaluationScheme
    {
    public:
        MLSEvaluationScheme()
        {
            setPrecision(Eigen::NumTraits<Scalar>::dummy_precision());
            setNIter(5);
        }

        MLSEvaluationScheme(unsigned int nIter) : MLSEvaluationScheme() { setNIter(nIter); }

        MLSEvaluationScheme(unsigned int nIter, Scalar eps) : MLSEvaluationScheme(nIter)
        {
            setPrecision(Eigen::NumTraits<Scalar>::dummy_precision());
        }

        /**
         * \brief Sets precision for early stopping
         *
         * \param _eps
         */
        void setPrecision(Scalar _eps) { eps = _eps; }

        /**
         * \brief Set maximum number of iterations
         *
         * \param _nIter Maximum number of iterations
         */
        void setNIter(unsigned int _nIter)
        {
            // At least 1, otherwise no computation would be performed
            nIter = nIter;
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam ComputeObject ComputeObject type
         * \tparam ItB Begin iterator type
         * \tparam ItE End iterator type
         * \tparam Project Projection functor type
         *
         * \param _co The fitting object
         * \param _beg The begining of point range
         * \param _end The end of point range
         * \param _p Projection functor
         *
         * \return The result of the fit
         */
        template <typename ComputeObject, typename ItB, typename ItE, typename Project = DirectProjectionOperator>
        PONCA_MULTIARCH inline FIT_RESULT compute(ComputeObject& _co, const ItB& _itb, const ItE& _ite,
                                                  const Project& _p = Project{}) const
        {
            return computeMLSImpl(_co, [&]() { return _co.compute(_itb, _ite); }, _p);
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam ComputeObject ComputeObject type
         * \tparam Container Container of points
         * \tparam Project Projection functor type
         *
         * \param _co The fitting object
         * \param _container The point container
         * \param _p Projection functor
         *
         * \return The result of the fit
         */
        template <typename ComputeObject, typename PointContainer, typename Project = DirectProjectionOperator>
        PONCA_MULTIARCH inline FIT_RESULT compute(ComputeObject& _co, const PointContainer& container,
                                                  const Project& _p = Project{})
        {
            return compute(_co, std::begin(container), std::end(container), _p);
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam ComputeObject ComputeObject type
         * \tparam IdxRange Range index type
         * \tparam PointContainer Point container (must provide random access)
         * \tparam Project Projection functor type
         *
         * \param _co The fitting object
         * \param _range The container of indices
         * \param _container The point container
         * \param _p Projection functor
         *
         * \return The result of the fit
         */
        template <typename ComputeObject, typename IdxRange, typename PointContainer,
                  typename Project = DirectProjectionOperator>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(ComputeObject& _co, const IdxRange& _range,
                                                         const PointContainer& _container,
                                                         const Project& _p = Project{}) const
        {
            return computeMLSImpl(_co, [&]() { return _co.computeWithIds(_range, _container); }, _p);
        }

    public:
        Scalar eps         = Eigen::NumTraits<Scalar>::dummy_precision();
        unsigned int nIter = 5;

    private:
        /*!
         * \brief Computes the fit using the MLS iteration process.
         *
         * The position of the projected point is outputted within getNeighborFilter().evalPos()
         *
         * \tparam ComputeObject The fitting type
         * \tparam Func The compute procedure
         * \tparam Project Projection functor type
         *
         * \param _co The fitting object
         * \param _compute The procedure to compute estimator
         * \param _p Projection functor
         *
         * \return The result of the fit
         */
        template <typename ComputeObject, typename Func, typename Project = DirectProjectionOperator>
        PONCA_MULTIARCH inline FIT_RESULT computeMLSImpl(ComputeObject& _co, Func&& _compute,
                                                         const Project& _p = Project{}) const
        {
            FIT_RESULT res = UNDEFINED;
            auto filter    = _co.getNeighborFilter();
            auto lastPos   = filter.evalPos();

            for (unsigned int mm = 0; mm < nIter; ++mm)
            {
                filter.changeNeighborhoodFrame(lastPos);
                _co.setNeighborFilter(filter);
                res = _compute();

                if (_co.isStable())
                {
                    auto newPos = _p(_co, lastPos);
                    if (newPos.isApprox(lastPos, eps))
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
