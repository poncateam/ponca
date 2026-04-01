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
     * \brief Computes the fit using the MLS iteration process
     *
     * The position of the projected point is outputted within getNeighborFilter().evalPos()
     *
     * \tparam Scalar scalar type
     */
    template <typename Scalar>
    class MLS
    {
    public:
        /**
         * \brief Default constructor
         *
         * - Precision is Eigen dummy_precision
         * - 5 iterations maximum
         */
        MLS()
        {
            setPrecision(Eigen::NumTraits<Scalar>::dummy_precision());
            setNIter(5);
        }

        /**
         * \brief Constructor from number of iterations
         *
         * \param _nIter Maximum number of iterations.
         */
        MLS(unsigned int _nIter) : MLS() { setNIter(_nIter); }

        /**
         * \brief Constructor from number of iterations and precision
         *
         * \param _nIter Maximum number of iterations.
         * \param _eps Precision for early stopping.
         */
        MLS(unsigned int _nIter, Scalar _eps) : MLS(_nIter)
        {
            setPrecision(Eigen::NumTraits<Scalar>::dummy_precision());
        }

        /**
         * \brief Sets precision for early stopping
         *
         * \param _eps
         */
        void setPrecision(Scalar _eps) { m_eps = _eps; }

        /**
         * \brief Set maximum number of iterations
         *
         * \param _nIter Maximum number of iterations
         */
        void setNIter(unsigned int _nIter)
        {
            // At least 1, otherwise no computation would be performed
            _nIter = std::max(m_nIter, 1u);
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam ItB Begin iterator type
         * \tparam ItE End iterator type
         *
         * \param _fit The fitting object
         * \param _beg The begining of point range
         * \param _end The end of point range
         *
         * \return The result of the fit
         */
        template <typename Fit, typename ItB, typename ItE>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& _fit, const ItB& _itb, const ItE& _ite) const
        {
            return computeMLSImpl(_fit, [&]() { return _fit.compute(_itb, _ite); });
        }

        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam Container Container of points
         *
         * \param _fit The fitting object
         * \param _container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& fit, const PointContainer& container) const
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
         * \param _fit The fitting object
         * \param _range The container of indices
         * \param _container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename IdxRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(Fit& _fit, const IdxRange& _range,
                                                         const PointContainer& _container) const
        {
            return computeMLSImpl(_fit, [&]() { return _fit.computeWithIds(_range, _container); });
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
         * \param _fit The fitting object
         * \param _compute The procedure to compute estimator
         *
         * \return The result of the fit
         */
        template <typename Fit, typename Func>
        PONCA_MULTIARCH inline FIT_RESULT computeMLSImpl(Fit& _fit, Func&& _compute) const
        {
            FIT_RESULT res = UNDEFINED;
            auto filter    = _fit.getNeighborFilter();
            auto lastPos   = filter.evalPos();

            for (unsigned int mm = 0; mm < m_nIter; ++mm)
            {
                filter.changeNeighborhoodFrame(lastPos);
                _fit.setNeighborFilter(filter);
                res = _compute();

                if (_fit.isStable())
                {
                    auto newPos = _fit.project(lastPos);
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

    /*!
     * \brief Computes the fit using the RIMLS iteration process \cite oztireli:2009:feature
     *
     * The position of the projected point is outputted within getNeighborFilter().evalPos()
     *
     * \tparam Scalar scalar type
     */
    template <typename Scalar, typename WeightKernel>
    class RIMLS
    {
    public:
        /**
         * \brief Default constructor
         *
         * - Default constructs weight kernel
         * - /!\ Uses sigma = 0.5. This should be scaled by the 'object scale'
         * - Precisions are set to Eigen dummy_precision
         * - Iterations are set to 5
         */
        RIMLS()
        {
            constexpr Scalar eps = Eigen::NumTraits<Scalar>::dummy_precision();

            setKernel(WeightKernel{});
            setSigma(0.5);

            setPrecisions(eps, eps);
            setIterations(5, 5);
        }

        /**
         * \see setSigma
         */
        RIMLS(Scalar _sigma) : RIMLS() { setSigma(_sigma); }

        /**
         * \brief Sets weight kernel
         *
         * This is named w in \cite oztireli:2009:feature. The code will use the
         * `f` method internally.
         *
         * \param _kernel The new weight kernel
         */
        void setKernel(WeightKernel _kernel) { m_W = std::move(_kernel); }

        /**
         * \brief Sets kernel scale factor
         *
         * The recommanded value is 0.5h, where h depends on the object scale (\cite oztireli:2009:feature)
         *
         * \param _sigma The kernel scale factor
         */
        void setSigma(Scalar _sigma) { m_sigma = _sigma; }

        /**
         * \brief Sets precision for fitting position and scalar field
         *
         * \param _epsX precision for the position
         * \param _epsF precision for the scalar field
         */
        void setPrecisions(Scalar _epsX, Scalar _epsF)
        {
            m_epsX = _epsX;
            m_epsF = _epsF;
        }

        /**
         * \brief Sets iterations for fitting position and scalar field
         *
         * /!\ The maximum complexity is _fIt * _xIt * #points
         *
         * \param _xIt iterations for the position
         * \param _fIt iterations for the scalar field
         */
        void setIterations(unsigned int _xIt, unsigned int _fIt)
        {
            m_maxIterX = _xIt;
            m_maxIterF = _fIt;
        }

    public:
        /**
         * \copydoc computeMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam ItB Begin iterator type
         * \tparam ItE End iterator type
         *
         * \param _fit The fitting object
         * \param _beg The begining of point range
         * \param _end The end of point range
         *
         * \return The result of the fit
         */
        template <typename Fit, typename ItB, typename ItE>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& _fit, const ItB& _itb, const ItE& _ite) const
        {
            return computeRIMLSImpl(
                _fit, _itb, _ite, [](const auto& _it) { return *_it; }, [&]() { return _fit.compute(_itb, _ite); });
        }

        /**
         * \copydoc computeRIMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam Container Container of points
         *
         * \param _fit The fitting object
         * \param _container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT compute(Fit& _fit, const PointContainer& _container) const
        {
            return compute(_fit, std::begin(_container), std::end(_container));
        }

        /**
         * \copydoc computeRIMLSImpl
         *
         * \tparam Fit Fit type
         * \tparam IdxRange Range index type
         * \tparam PointContainer Point container (must provide random access)
         *
         * \param _fit The fitting object
         * \param _range The container of indices
         * \param _container The point container
         *
         * \return The result of the fit
         */
        template <typename Fit, typename IdxRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(Fit& _fit, const IdxRange& _range,
                                                         const PointContainer& _container) const
        {
            // TODO: when availaible use cuda::std::views instead of a lambda
            return computeRIMLSImpl(
                _fit, std::begin(_range), std::end(_range), [&](const auto& _it) { return _container[*_it]; },
                [&]() { return _fit.computeWithIds(_range, _container); });
        }

    protected:
        WeightKernel m_W; //< Weight kernel

        Scalar m_sigma; //< Kernel scale factor

        Scalar m_epsX;           //< Precision for fitting the position
        unsigned int m_maxIterX; //< Max number of iteration for fitting the position

        Scalar m_epsF;           //< Precision for fitting the scalar field
        unsigned int m_maxIterF; //< Max number of iteration for fitting the scalar field

        /**
         * \brief Computes the fit using the RIMLS scheme
         *
         * The method was published in \cite Oztireli:2009
         *
         * The maximum total number of iteration is _maxIterF * _maxIterX. It is advised to keep
         * both values as low as possible in order to avoid high computationnal times.
         *
         * \tparam Fit Fit type
         * \tparam ItB Iterator type
         * \tparam ItE Iterator type
         * \tparam Func Iterator-to-point conversion function
         *
         * \param _fit Fit object
         * \param _begin Begining of range
         * \param _end End of range
         * \param _extractPoint Iterator-to-point converter
         */
        template <typename Fit, typename ItB, typename ItE, typename Func, typename Compute>
        PONCA_MULTIARCH FIT_RESULT computeRIMLSImpl(Fit& _fit, ItB _beg, ItE _end, Func&& _extractPoint,
                                                    Compute&& _compute) const
        {
            using VectorType = typename Fit::VectorType;

            FIT_RESULT res = _compute();

            if (res != STABLE)
                return res;

            PONCA_MULTIARCH_STD_MATH(abs);
            PONCA_MULTIARCH_STD_MATH(sqrt);

            Scalar sigma = m_sigma;

            auto filter = _fit.getNeighborFilter();

            Scalar f         = _fit.potential();
            VectorType x     = filter.evalPos();
            VectorType gradF = _fit.primitiveGradient();

            unsigned int xIt = 0;
            unsigned int fIt = 0;

            Scalar xRes = Eigen::NumTraits<Scalar>::infinity();
            Scalar fRes = Eigen::NumTraits<Scalar>::infinity();

            while (xIt < m_maxIterX && xRes < m_epsX)
            {
                while (fIt < m_maxIterF && fRes < m_epsF)
                {
                    Scalar sumF      = 0;
                    Scalar sumW      = 0;
                    VectorType sumGW = VectorType::Zero();
                    VectorType sumGF = VectorType::Zero();
                    VectorType sumN  = VectorType::Zero();

                    for (auto it = _beg; it != _end; ++it)
                    {
                        const auto pt       = _extractPoint(it);
                        const VectorType px = x - pt.pos();
                        const Scalar fx     = px.dot(pt.normal());

                        Scalar alpha = 1;
                        if (xIt > 0)
                            alpha = m_W.f((fx - f) / sigma) * m_W.f((pt.normal() - gradF).norm() / sigma);

                        const Scalar w      = alpha * filter(px).first;
                        const VectorType gw = 2 * alpha * px * filter.scaledw(x, pt.pos());

                        sumW += w;
                        sumGW += gw;
                        sumF += w * fx;
                        sumGF += gw * fx;
                        sumN += w * pt.normal();
                    }

                    // [Oztireli 2009] defines another residual formula which depends
                    // on all neighbors and requires a storage of twice this size.
                    // Here, we replaced this with an absolute convergence criterion on
                    // the scalar field values for faster and easier computation.
                    const Scalar newF = sumF / sumW;
                    fRes              = abs(newF - f);

                    f     = newF;
                    gradF = (sumGF - f * sumGW + sumN) / sumW;

                    ++fIt;
                }

                const VectorType step = f * gradF;
                xRes                  = step.squaredNorm();

                x -= step;
                ++xIt;
            }

            filter.changeNeighborhoodFrame(x);

            return STABLE;
        }
    };
} // namespace Ponca
