/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "../concepts.h"
#include "../defines.h"
#include "../enums.h"

#include "./project.h"

#include "Eigen/Core"

namespace Ponca
{
    /*!
     * \brief Computes the fit using the Moving Least Squares approach.
     * The projection operator can be customized, see \cite Alexa:2004:projection.
     *
     * The position of the projected point is outputted within getNeighborFrame().center()
     *
     * \tparam Scalar scalar type
     */
    template <typename Scalar>
    struct MLSEvaluationScheme
    {
        MLSEvaluationScheme(unsigned int _nIter = nIterDefault, Scalar _eps = epsDefault) : nIter(_nIter), eps(_eps) {}

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

        /// \brief Default epsilon value for stopping MLS iterations
        static constexpr Scalar epsDefault = Eigen::NumTraits<Scalar>::dummy_precision();
        /// \brief Default maximum number of MLS iterations
        static constexpr int nIterDefault = 5;
        /// \brief Epsilon value for stopping MLS iterations
        Scalar eps = epsDefault;
        /// \brief Maximum number of MLS iterations
        unsigned int nIter = nIterDefault;

    private:
        /*!
         * \brief Computes the fit using the MLS iteration process.
         *
         * The position of the projected point is outputted within getNeighborFrame().center()
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
            requires ProvidesBasketUnitBase<ComputeObject> && ProvidesImplicitPrimitive<ComputeObject>
        PONCA_MULTIARCH inline FIT_RESULT computeMLSImpl(ComputeObject& _co, Func&& _compute,
                                                         const Project& _p = Project{}) const
        {
            FIT_RESULT res = UNDEFINED;
            auto& frame    = _co.getNeighborFrame();
            auto lastPos   = frame.center();

            for (unsigned int mm = 0; mm < nIter; ++mm)
            {
                frame.changeNeighborhoodFrame(lastPos);
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
