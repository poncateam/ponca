/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../Common/Assert.h"
#include "./defines.h"
#include "./frame.h"
#include PONCA_MULTIARCH_INCLUDE_CU_STD(utility)

namespace Ponca
{
    /*!
        \brief Weight neighbors according to the Euclidean distance between a query and a reference position

        The evaluation position is set in the constructor. All the queries are expressed in global system, and converted
        internally to relatively to the evaluation position using #CenteredNeighborhoodFrame.

        \tparam DataPoint Type of input points.
        \tparam WeightKernel 1d function used to compute the weight depending on the distance between query point and
       the basis center. If `WeightKernel::isCompact == true`, the distance to the basis center is checked against the
       scale parameter, and any point whose distance is larger than the scale will be assigned with a weight of \f$0\f$.
       For non-compact (ie. global) kernels, the weights are computed for all points.

        \see operator()

        \warning DistWeightFunc assumes that the evaluation scale t is strictly positive, but the valus is not checked
    */
    template <class DataPoint, class WeightKernel>
    class DistWeightFilter : public CenteredNeighborhoodFrame<DataPoint>
    {
    public:
        /*! \brief Scalar type from DataPoint */
        using Scalar = typename DataPoint::Scalar;
        /*! \brief Vector type from DataPoint */
        using VectorType = typename DataPoint::VectorType;
        /*! \brief Matrix type from DataPoint */
        using MatrixType = typename DataPoint::MatrixType;
        /*! \brief Return type of the method #w() */
        using WeightReturnType = PONCA_MULTIARCH_CU_STD_NAMESPACE(pair)<Scalar, VectorType>;
        /// \brief Frame used to express the neighbors locally
        using NeighborhoodFrame = CenteredNeighborhoodFrame<DataPoint>;
        /// \brief Flag indicating if the weighting kernel is compact of not
        static constexpr bool isCompact = WeightKernel::isCompact;

        /*!
            \brief Constructor that defines the current evaluation scale
            \warning t > 0
        */
        PONCA_MULTIARCH inline DistWeightFilter(const VectorType& _evalPos = VectorType::Zero(),
                                                const Scalar& _t           = Scalar(1.))
            : NeighborhoodFrame(_evalPos), m_t(_t)
        {
            PONCA_ASSERT(_t > Scalar(0));
        }

        ///! \copydoc DistWeightFunc
        PONCA_MULTIARCH inline DistWeightFilter(const DataPoint& _evalPoint, const Scalar& _t = Scalar(1.))
            : NeighborhoodFrame(_evalPoint.pos()), m_t(_t)
        {
            PONCA_ASSERT(_t > Scalar(0));
        }

        /*!
            \brief Compute the weight of the given query with respect to its coordinates.

            \param _q Query in global coordinate system

            As the query \f$\mathbf{q}\f$ is expressed in global coordinate, it is
            first converted to the centered basis. Then, the WeightKernel is directly
            applied to the norm of its coordinates with respect to the current scale  \f$ t \f$ :

            \f$ w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) \f$

            \see convertToLocalBasis
            \return The computed weight + the point expressed in local basis
        */
        PONCA_MULTIARCH inline WeightReturnType operator()(const DataPoint& q) const;

        /*!
            \brief First order derivative in space (for each spatial dimension \f$\mathsf{x})\f$

            \param _q Query in global coordinate system

            \f$ \frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}}
            \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t})
            = \frac{\mathbf{q}}{t\left|q\right|}  \nabla{w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t})}  \f$

            where \f$ \left|\mathbf{q}_\mathsf{x}\right| \f$ represents the norm of the
            query coordinates expressed in centered basis,
            for each spatial dimensions \f$ \mathsf{x}\f$.

            \warning Requires \f$\nabla w(x)\f$ to be valid
        */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType spacedw(const VectorType& _q,
                                                                const DataPoint& /*attributes*/) const;

        /*!
            \brief Second order derivative in space (for each spatial dimension \f$\mathsf{x})\f$

            \param _q Query in global coordinate system

            \f$ \frac{\delta^2 \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}^2}
            \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
            \left(\frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}}\right)^2
            \nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) =
            \frac{1}{t\left|\mathbf{q}_\mathsf{x}\right|} \left( I_d -
           \frac{\mathbf{q}_\mathsf{x}\mathbf{q}_\mathsf{x}^T}{\left|\mathbf{q}_\mathsf{x}\right|^2}\right)
            \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
            \frac{\mathbf{q}_\mathsf{x}\mathbf{q}_\mathsf{x}^T}{t^2\left|\mathbf{q}_\mathsf{x}\right|^2}
            \nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) \f$

            where \f$ \left|\mathbf{q}_\mathsf{x}\right| \f$ represents the norm of the
            query coordinates expressed in centered basis,
            for each spatial dimensions \f$ \mathsf{x}\f$.

            \warning Requires \f$\nabla^2 w(x)\f$ to be valid
        */
        PONCA_MULTIARCH [[nodiscard]] inline MatrixType spaced2w(const VectorType& _q,
                                                                 const DataPoint& /*attributes*/) const;

        /*!
            \brief First order derivative in scale  \f$t\f$

            \param _q Query in global coordinate system

            \f$ \frac{\delta \frac{\left|\mathbf{q}\right|}{t}}{\delta t}
            \nabla w(\frac{\left|\mathbf{q}\right|}{t})
            = - \frac{\left|\mathbf{q}\right|}{t^2} \nabla{w(\frac{\left|\mathbf{q}\right|}{t})} \f$

            where \f$ \left|\mathbf{q}\right| \f$ represents the norm of the
            query coordinates expressed in centered basis.

            \warning Requires \f$\nabla w(x)\f$ to be valid
        */
        PONCA_MULTIARCH [[nodiscard]] inline Scalar scaledw(const VectorType& _q,
                                                            const DataPoint& /*attributes*/) const;

        /*!
            \brief Second order derivative in scale  \f$t\f$

            \param _q Query in global coordinate system

            \f$ \frac{\delta^2 \frac{\left|\mathbf{q}\right|}{t}}{\delta t^2}
            \nabla w(\frac{\left|\mathbf{q}\right|}{t}) +
            \left(\frac{\delta \frac{\left|\mathbf{q}\right|}{t}}{\delta t}\right)^2
            \nabla^2 w(\frac{\left|\mathbf{q}\right|}{t}) =
            \frac{2\left|\mathbf{q}\right|}{t^3} \nabla{w(\frac{\left|\mathbf{q}\right|}{t})} +
            \frac{\left|\mathbf{q}\right|^2}{t^4} \nabla^2{w(\frac{\left|\mathbf{q}\right|}{t})}
            \f$

            where \f$ \left|\mathbf{q}\right| \f$ represents the norm of the
            query coordinates expressed in centered basis.

            \warning Requires \f$\nabla^2 w(x)\f$ to be valid
        */
        PONCA_MULTIARCH [[nodiscard]] inline Scalar scaled2w(const VectorType& _q,
                                                             const DataPoint& /*attributes*/) const;

        /*!
            \brief Cross derivative in scale \f$t\f$ and in space (for each spatial dimension \f$\mathsf{x})\f$

            \param _q Query in global coordinate system

            \f$ \frac{\delta^2 \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta t\ \delta \mathsf{x}}
            \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
            \frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}}
            \frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta t}
            \nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) =
            -\frac{\mathbf{q}_\mathsf{x}}{t^2}
            \left( \frac{1}{\left|\mathbf{q}_\mathsf{x}\right|}\nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
            \frac{1}{t}\nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) \right)\f$

            where \f$ \left|\mathbf{q}_\mathsf{x}\right| \f$ represents the norm of the
            query coordinates expressed in centered basis.

            \warning Requires \f$\nabla^2 w(x)\f$ to be valid
        */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType scaleSpaced2w(const VectorType& _q,
                                                                      const DataPoint& /*attributes*/) const;

        /*! \brief Access to the evaluation scale set during the initialization */
        PONCA_MULTIARCH [[nodiscard]] inline Scalar evalScale() const { return m_t; }

        PONCA_MULTIARCH [[nodiscard]] inline NeighborhoodFrame& frame() { return *this; }
        PONCA_MULTIARCH [[nodiscard]] inline const NeighborhoodFrame& frame() const { return *this; }

    protected:
        Scalar m_t;        /*!< \brief Evaluation scale */
        WeightKernel m_wk; /*!< \brief 1D function applied to weight queries */

    }; // class DistWeightFunc

    namespace internal
    {
        /*!
            \brief Weighting function that set uniform weight to all samples

            In contrast to DistWeightFunc with ConstantWeight, it does not check for scale range.
            \tparam _NeighborhoodFrame Base NeighborhoodFrame used to performs (or not) local basis conversion and
           maintain computation accuracy
        */
        template <class DataPoint, template <typename> typename _NeighborhoodFrame>
        class NoWeightFilterBase : public _NeighborhoodFrame<DataPoint>
        {
        public:
            /*! \brief Scalar type from DataPoint */
            using Scalar = typename DataPoint::Scalar;
            /*! \brief Vector type from DataPoint */
            using VectorType = typename DataPoint::VectorType;
            /*! \brief Matrix type from DataPoint */
            using MatrixType = typename DataPoint::MatrixType;
            /*! \brief Return type of the method #w() */
            using WeightReturnType = PONCA_MULTIARCH_CU_STD_NAMESPACE(pair)<Scalar, VectorType>;

            using NeighborhoodFrame = _NeighborhoodFrame<DataPoint>;

            /*!
                \brief Default constructor. All parameters are ignored (kept for API compatibility with DistWeightFunc.
            */
            PONCA_MULTIARCH inline NoWeightFilterBase(const VectorType& v = VectorType::Zero(), Scalar = 0)
                : NeighborhoodFrame(v)
            {
            }

            ///! \copydoc NoWeightFuncBase
            PONCA_MULTIARCH inline NoWeightFilterBase(const DataPoint& v, Scalar = 0) : NoWeightFilterBase(v.pos()) {}

            /*!
                \brief Compute the weight of the given query, which is always $1$
                \param _q Query in global coordinate system
            */
            PONCA_MULTIARCH inline WeightReturnType operator()(const DataPoint& _q) const
            {
                return {Scalar(1), NeighborhoodFrame::convertToLocalBasis(_q.pos())};
            }

            /*!
                \brief First order derivative in space (for each spatial dimension \f$\mathsf{x})\f$, which are always
               $0$
                \param _q Query in global coordinate system
            */
            PONCA_MULTIARCH [[nodiscard]] inline VectorType spacedw(const VectorType& /*_q*/,
                                                                    const DataPoint& /*attributes*/) const
            {
                return VectorType::Zeros();
            }

            /*!
                \brief Second order derivative in space (for each spatial dimension \f$\mathsf{x})\f$, which are always
               $0$
                \param _q Query in global coordinate
            */
            PONCA_MULTIARCH [[nodiscard]] inline MatrixType spaced2w(const VectorType& /*_q*/,
                                                                     const DataPoint& /*attributes*/) const
            {
                return MatrixType::Zeros();
            }

            /*!
                \brief First order derivative in scale  \f$t\f$, which are always $0$
                \param _q Query in global coordinate
            */
            PONCA_MULTIARCH [[nodiscard]] inline Scalar scaledw(const VectorType& /*_q*/,
                                                                const DataPoint& /*attributes*/) const
            {
                return Scalar(0);
            }

            /*!
                \brief Second order derivative in scale  \f$t\f$, which are always $0$
                \param _q Query in global coordinate
            */
            PONCA_MULTIARCH [[nodiscard]] inline Scalar scaled2w(const VectorType& /*_q*/,
                                                                 const DataPoint& /*attributes*/) const
            {
                return Scalar(0);
            }

            /*!
                \brief Cross derivative in scale \f$t\f$ and in space (for each spatial dimension \f$\mathsf{x})\f$,
               which are always $0$
                \param _q Query in global coordinate
            */
            PONCA_MULTIARCH [[nodiscard]] inline VectorType scaleSpaced2w(const VectorType& /*_q*/,
                                                                          const DataPoint& /*attributes*/) const
            {
                return VectorType::Zeros();
            }

            PONCA_MULTIARCH [[nodiscard]] inline NeighborhoodFrame& frame() { return *this; }
            PONCA_MULTIARCH [[nodiscard]] inline const NeighborhoodFrame& frame() const { return *this; }
        }; // class NoWeightFuncBase

        template <typename DataPoint>
        inline constexpr typename DataPoint::VectorType ConvertDataPointToNormal(const DataPoint& pt)
        {
            return pt.normal();
        }
    } // namespace internal

    /*!
     * This class extends a NeighborFilter class to also store additional data, for use outside the
     * scope of this class.
     *
     * \tparam DataType The type of data to store along the filter
     * \tparam NeighborFilter Any NeighborFilter type (e.g., NoWeightFilter or DistWeightFilter<ConstantWeightKernel>)
     * \tparam DataConverter A function that construct a DataType from the DataPoint
     */
    template <class DataPoint, typename DataType, typename NeighborFilter, auto DataConverter>
    class FilterWithAttributes : public NeighborFilter
    {
    public:
        using Base = NeighborFilter;
        /*! \brief Scalar type from DataPoint */
        using Scalar = typename DataPoint::Scalar;
        /*! \brief Vector type from DataPoint */
        using VectorType = typename DataPoint::VectorType;
        /*! \brief Matrix type from DataPoint */
        using MatrixType = typename DataPoint::MatrixType;
        /*! \brief Return type of the method #w() */
        using WeightReturnType = PONCA_MULTIARCH_CU_STD_NAMESPACE(pair)<Scalar, VectorType>;

        /*!
            \brief Constructor that defines the current evaluation scale
            \warning t > 0
        */
        PONCA_MULTIARCH inline FilterWithAttributes(const VectorType& _evalPos = VectorType::Zero(),
                                                    const Scalar& _t = Scalar(1), const DataType& _data = DataType{})
            : Base(_evalPos, _t), m_data(_data)
        {
        }

        /*!
            \brief Constructor that defines the current evaluation scale and construct the DataType from the _evalPoint
            using the DataConverter
            \warning t > 0
        */
        PONCA_MULTIARCH inline FilterWithAttributes(const DataPoint& _evalPoint, const Scalar& _t = Scalar(0))
            : Base(_evalPoint.pos(), _t), m_data(DataConverter(_evalPoint))
        {
        }

        /*! \brief Access to the evaluation normal set during the initialization */
        PONCA_MULTIARCH inline const DataType& data() const { return m_data; }

    protected:
        DataType m_data; /*!< \brief Evaluation normal */
    }; // class FilterWithAttributes

    template <class DataPoint>
    /// \brief Weighting function that set uniform weight to all samples, but transform neighbors coordinates to local
    /// frame
    /// \see internal::NoWeightFuncBase
    struct NoWeightFilter : public internal::NoWeightFilterBase<DataPoint, CenteredNeighborhoodFrame>
    {
        using internal::NoWeightFilterBase<DataPoint, CenteredNeighborhoodFrame>::NoWeightFilterBase;
    };

    template <class DataPoint>
    /// \brief Weighting function that set uniform weight to all samples and keep neighbors coordinates in global frame
    /// \see internal::NoWeightFuncBase
    struct NoWeightFilterGlobal : public internal::NoWeightFilterBase<DataPoint, GlobalNeighborhoodFrame>
    {
        using internal::NoWeightFilterBase<DataPoint, GlobalNeighborhoodFrame>::NoWeightFilterBase;
    };

#include "weightFilter.hpp"

} // namespace Ponca
