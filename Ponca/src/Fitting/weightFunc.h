/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include PONCA_MULTIARCH_INCLUDE_CU_STD(utility)

namespace Ponca
{
/*!
 * \brief NeighborhoodFrame that express 3d points relatively to a prescribed center.
 *
 * This class is useful to get all coordinates centered around a point, which ease weights computation, and limits
 * issues with big numbers and rounding errors.
 *
 * Express points \f$\mathbf{x}\f$ relatively to a center \f$\mathbf{p}\f$, ie.
 * \f$\mathbf{x}'=\mathbf{x}-\mathbf{p}\f$.
 * This frame does not apply rotation.
 *
 * @tparam DataPoint Point type used for computation
 */
template<class DataPoint>
class CenteredNeighborhoodFrame {
public:
    /*! \brief Scalar type from DataPoint */
    using Scalar =  typename DataPoint::Scalar;
    /*! \brief Vector type from DataPoint */
    using VectorType =  typename DataPoint::VectorType;

    /// \brief Flag indicating that this class modifies the coordinates when passing from global to local
    static constexpr bool hasLocalFrame = true;

    PONCA_MULTIARCH inline explicit CenteredNeighborhoodFrame(const VectorType & _evalPos = VectorType::Zero())
    : m_p(_evalPos) {}

    /// \brief Change neighborhood frame (move basis center)
    PONCA_MULTIARCH inline void changeNeighborhoodFrame(const VectorType& _newEvalPos) { m_p = _newEvalPos; };

    /*!
     * \brief Convert query from local to global coordinate system, such as \f$\mathbf{x}=\mathbf{x}'+\mathbf{p}\f$.
     * @param _q Position expressed relatively to the basis center
     * @return Position in the global coordinate system
     */
    PONCA_MULTIARCH inline VectorType convertToGlobalBasis(const VectorType& _q) const { return _q + m_p; }

    /*!
     * \brief Convert query from global to local coordinate system, such as \f$\mathbf{x}'=\mathbf{x}-\mathbf{p}\f$.
     *
     * @param _q Input point in global coordinate system
     * @return Position expressed relatively to the basis center
     */
    PONCA_MULTIARCH inline VectorType convertToLocalBasis(const VectorType& _q) const { return _q - m_p; }

    /*!
     * \brief Get access to the stored points of evaluation
     * @return Position of the local basis center
     */
    PONCA_MULTIARCH inline VectorType evalPos() const { return m_p; }

private:
    VectorType m_p;  /*!< \brief basis center */
};
/*!
 * \brief NeighborhoodFrame that keep points in the global frame without applying any transformation
 *
 * This class is useful to compute
 * Express points \f$\mathbf{x}\f$ relatively to a center \f$\mathbf{p}\f$, ie.
 * \f$\mathbf{x}'=\mathbf{x}-\mathbf{p}\f$.
 * This frame does not apply rotation.
 *
 * @tparam DataPoint Point type used for computation
 */
template<class DataPoint>
class GlobalNeighborhoodFrame {
public:
    /*! \brief Scalar type from DataPoint */
    using Scalar =  typename DataPoint::Scalar;
    /*! \brief Vector type from DataPoint */
    using VectorType =  typename DataPoint::VectorType;

    /// \brief Flag indicating that this class does not modify the coordinates when passing from global to local
    static constexpr bool hasLocalFrame = false;

    PONCA_MULTIARCH inline explicit GlobalNeighborhoodFrame(const VectorType & /*_evalPos*/ = VectorType::Zero()) {}

    /// \brief Change neighborhood frame (has no effect for global basis)
    PONCA_MULTIARCH inline void changeNeighborhoodFrame(const VectorType& /*_newEvalPos*/) {};

    /*!
     * \brief Convert position from local to global coordinate system : does nothing as this is global frame
     * @param _q Position in local coordinate
     * @return _q
     */
    PONCA_MULTIARCH inline const VectorType& convertToGlobalBasis(const VectorType& _q) const { return _q; }

    /*!
     * \brief Convert query from global to local coordinate system : does nothing as this is global frame
     * @param _q Query in global coordinate
     * @return _q
     */
    PONCA_MULTIARCH inline const VectorType& convertToLocalBasis(const VectorType& _q) const { return _q; }
};


/*!
    \brief Weight neighbors according to the euclidean distance between a query and a reference position

    The evaluation position is set in the constructor. All the queries are expressed in global system, and converted
    internally to relatively to the evaluation position using #CenteredNeighborhoodFrame.

    \tparam DataPoint Type of input points.
    \tparam WeightKernel 1d function used to compute the weight depending on the distance between query point and the
    basis center. If `WeightKernel::isCompact == true`, the distance to the basis center is checked against the scale
    parameter, and any point whose distance is larger than the scale will be assigned with a weight of \f$0\f$. For
    non-compact (ie. global) kernels, the weights are computed for all points.

    \see operator()

    \warning DistWeightFunc assumes that the evaluation scale t is strictly positive, but the valus is not checked
*/
template <class DataPoint, class WeightKernel>
class DistWeightFunc : public CenteredNeighborhoodFrame<DataPoint>
{
public:
    /*! \brief Scalar type from DataPoint */
    using Scalar =  typename DataPoint::Scalar;
    /*! \brief Vector type from DataPoint */
    using VectorType =  typename DataPoint::VectorType;
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
    PONCA_MULTIARCH inline DistWeightFunc(const VectorType & _evalPos = VectorType::Zero(),
                                              const Scalar& _t = Scalar(1.))
    : NeighborhoodFrame(_evalPos), m_t(_t)
    {
        //\todo manage that assert on __host__ and __device__
        //assert(_t > Scalar(0));
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
    PONCA_MULTIARCH inline VectorType spacedw(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;


    /*!
        \brief Second order derivative in space (for each spatial dimension \f$\mathsf{x})\f$

        \param _q Query in global coordinate system

        \f$ \frac{\delta^2 \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}^2}
        \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
        \left(\frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}}\right)^2
        \nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) =
        \frac{1}{t\left|\mathbf{q}_\mathsf{x}\right|} \left( I_d - \frac{\mathbf{q}_\mathsf{x}\mathbf{q}_\mathsf{x}^T}{\left|\mathbf{q}_\mathsf{x}\right|^2}\right)
        \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) +
        \frac{\mathbf{q}_\mathsf{x}\mathbf{q}_\mathsf{x}^T}{t^2\left|\mathbf{q}_\mathsf{x}\right|^2}
        \nabla^2 w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) \f$

        where \f$ \left|\mathbf{q}_\mathsf{x}\right| \f$ represents the norm of the
        query coordinates expressed in centered basis,
        for each spatial dimensions \f$ \mathsf{x}\f$.

        \warning Requires \f$\nabla^2 w(x)\f$ to be valid
    */
    PONCA_MULTIARCH inline MatrixType spaced2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

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
    PONCA_MULTIARCH inline Scalar scaledw(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

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
    PONCA_MULTIARCH inline Scalar scaled2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

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
    PONCA_MULTIARCH inline VectorType scaleSpaced2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

    /*! \brief Access to the evaluation scale set during the initialization */
    PONCA_MULTIARCH inline Scalar evalScale() const { return m_t; }

protected:
    Scalar       m_t;  /*!< \brief Evaluation scale */
    WeightKernel m_wk; /*!< \brief 1D function applied to weight queries */

};// class DistWeightFunc


/*!
    \brief Base Weighting function that set uniform weight to all samples

    In contrast to DistWeightFunc with ConstantWeight, it does not check for scale range.
    \tparam _NeighborhoodFrame Base NeighborhoodFrame used to performs (or not) local basis conversion and maintain
    computation accuracy
*/
template <class DataPoint, template <typename>typename _NeighborhoodFrame>
class NoWeightFuncBase : public _NeighborhoodFrame<DataPoint>
{
public:
    /*! \brief Scalar type from DataPoint */
    using Scalar =  typename DataPoint::Scalar;
    /*! \brief Vector type from DataPoint */
    using VectorType =  typename DataPoint::VectorType;
    /*! \brief Matrix type from DataPoint */
    using MatrixType = typename DataPoint::MatrixType;
    /*! \brief Return type of the method #w() */
    using WeightReturnType = PONCA_MULTIARCH_CU_STD_NAMESPACE(pair)<Scalar, VectorType>;

    using NeighborhoodFrame = _NeighborhoodFrame<DataPoint>;

    /*!
        \brief Default constructor. All parameters are ignored (kept for API compatibility with DistWeightFunc.
    */
    PONCA_MULTIARCH inline NoWeightFuncBase(const VectorType& v = VectorType::Zero(), Scalar = 0)
    : NeighborhoodFrame(v){ }

    /*!
        \brief Compute the weight of the given query, which is always $1$
        \param _q Query in global coordinate system
    */
    PONCA_MULTIARCH inline WeightReturnType operator()(const DataPoint& _q) const
    {
        return {Scalar(1), NeighborhoodFrame::convertToLocalBasis(_q.pos())};
    }


    /*!
        \brief First order derivative in space (for each spatial dimension \f$\mathsf{x})\f$, which are always $0$
        \param _q Query in global coordinate system
    */
    PONCA_MULTIARCH inline VectorType spacedw(const VectorType& /*_q*/,
                                              const DataPoint&  /*attributes*/) const
    { return VectorType::Zeros(); }


    /*!
        \brief Second order derivative in space (for each spatial dimension \f$\mathsf{x})\f$, which are always $0$
        \param _q Query in global coordinate
    */
    PONCA_MULTIARCH inline MatrixType spaced2w(const VectorType& /*_q*/,
                                               const DataPoint&  /*attributes*/) const
    { return MatrixType::Zeros(); }

    /*!
        \brief First order derivative in scale  \f$t\f$, which are always $0$
        \param _q Query in global coordinate
    */
    PONCA_MULTIARCH inline Scalar scaledw(const VectorType& /*_q*/,
                                          const DataPoint&  /*attributes*/) const
    { return Scalar(0); }

    /*!
        \brief Second order derivative in scale  \f$t\f$, which are always $0$
        \param _q Query in global coordinate
    */
    PONCA_MULTIARCH inline Scalar scaled2w(const VectorType& /*_q*/,
                                           const DataPoint&  /*attributes*/) const
    { return Scalar(0); }

    /*!
        \brief Cross derivative in scale \f$t\f$ and in space (for each spatial dimension \f$\mathsf{x})\f$, which are
        always $0$
        \param _q Query in global coordinate
    */
    PONCA_MULTIARCH inline VectorType scaleSpaced2w(const VectorType& /*_q*/,
                                                    const DataPoint&  /*attributes*/) const
    { return VectorType::Zeros(); }
};// class DistWeightFuncBase

template <class DataPoint>
/// \brief
using NoWeightFunc = NoWeightFuncBase<DataPoint, CenteredNeighborhoodFrame>;

template <class DataPoint>
/// \brief
using NoWeightFuncGlobal = NoWeightFuncBase<DataPoint, GlobalNeighborhoodFrame>;
#include "weightFunc.hpp"

}// namespace Ponca
