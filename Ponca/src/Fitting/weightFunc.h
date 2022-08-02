/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

namespace Ponca
{
    /*!
     \brief Applies a Concept::WeightKernelConcept to a Concept::PointConcept query.
     The input query is expressed in global coordinate system.
     A local frame (ie. a position) can be set a initialization time, and the query will be automatically converted
     from global to local coordinate at evaluation.
   */
    template <class DataPoint, class WeightKernel>
    class LocalWeightFunc {

    public:
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;

        /*!
         * \brief Initialization method, called by the fitting procedure
         * @param _evalPos Basis center
         */
        PONCA_MULTIARCH inline void init( const VectorType& _evalPos ) {}

        /*! \brief Convert query from global to local coordinate system */
        PONCA_MULTIARCH inline VectorType convertToLocalBasis(const VectorType& _q) const;

        /*! \brief Apply the weight function to a query. */
        PONCA_MULTIARCH inline Scalar w(const VectorType& globalQuery,
                                        const DataPoint&  attributes) const {}

        /*! \brief Apply the weight function derived in space to a query. */
        PONCA_MULTIARCH inline VectorType spacedw(const VectorType& globalQuery,
                                                  const DataPoint&  attributes) const {}

        /*! \brief Apply the weight function derived in scale to a query. */
        PONCA_MULTIARCH inline Scalar scaledw(const VectorType& globalQuery,
                                              const DataPoint&  attributes) const {}

        /*! \brief Read access to the evaluation scale */
        PONCA_MULTIARCH inline Scalar evalScale() const {}
    };// class WeightFuncConcept

/*!
    \brief Weighting function based on the euclidean distance between a query and a reference position

    The evaluation position is set using the #init method. All the queries are expressed in global system, and the
    weighting function convert them to local coordinates (ie. relatively to the evaluation position).

    This class inherits BaseWeightFunc. It can be specialized for any DataPoint,
    and uses a generic 1D BaseWeightKernel.

    \inherit Concept::WeightFuncConcept

    \warning it assumes that the evaluation scale t is strictly positive

    \ingroup fitting
*/
template <class DataPoint, class WeightKernel>
class DistWeightFunc
{
public:
    /*! \brief Scalar type from DataPoint */
    using Scalar =  typename DataPoint::Scalar;
    /*! \brief Vector type from DataPoint */
    using VectorType =  typename DataPoint::VectorType;
    /*! \brief Matrix type from DataPoint */
    using MatrixType = typename DataPoint::MatrixType;
    /*! \brief Return type of the method #w() */
    using WeightReturnType = std::pair<Scalar, VectorType>;

    /*!
        \brief Constructor that defines the current evaluation scale
        \warning t > 0
    */
    PONCA_MULTIARCH inline DistWeightFunc(const Scalar& _t = Scalar(1.))
    : m_p(VectorType::Zero())
    {
        //\todo manage that assrt on __host__ and __device__
        //assert(_t > Scalar(0));
        m_t = _t;
    }

    /*!
     * \brief Initialization method, called by the fitting procedure
     * @param _evalPos Basis center
     */
    PONCA_MULTIARCH inline void init( const VectorType& _evalPos )
    {
        m_p = _evalPos;
    }

    PONCA_MULTIARCH inline const VectorType& basisCenter() const { return m_p; }

    /*!
     * \brief Convert query from global to local coordinate system
     * @param _q Query in global coordinate
     * @return Query expressed relatively to the basis center
     */
    PONCA_MULTIARCH inline VectorType convertToLocalBasis(const VectorType& _q) const;

    /*!
        \brief Compute the weight of the given query with respect to its coordinates.

        As the query \f$\mathbf{q}\f$ is expressed in global coordinate, it is
        first converted to the centered basis. Then, the WeightKernel is directly
        applied to the norm of its coordinates with respect to the current scale  \f$ t \f$ :

        \f$ w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}) \f$

        \see convertToLocalBasis
        \return The computed weight + the point expressed in local basis
    */
    PONCA_MULTIARCH inline WeightReturnType w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;


    /*!
        \brief First order derivative in space (for each spatial dimension \f$\mathsf{x})\f$

        \f$ \frac{\delta \frac{\left|\mathbf{q}_\mathsf{x}\right|}{t}}{\delta \mathsf{x}}
        \nabla w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t})
        = \frac{\mathbf{q}}{t\left|q\right|}  \nabla{w(\frac{\left|\mathbf{q}_\mathsf{x}\right|}{t})}  \f$

        where \f$ \left|\mathbf{q}_\mathsf{x}\right| \f$ represents the norm of the
        query coordinates expressed in centered basis,
        for each spatial dimensions \f$ \mathsf{x}\f$.
    */
    PONCA_MULTIARCH inline VectorType spacedw(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;


    /*!
        \brief Second order derivative in space (for each spatial dimension \f$\mathsf{x})\f$

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
    */
    PONCA_MULTIARCH inline MatrixType spaced2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

    /*!
        \brief First order derivative in scale  \f$t\f$

        \f$ \frac{\delta \frac{\left|\mathbf{q}\right|}{t}}{\delta t}
        \nabla w(\frac{\left|\mathbf{q}\right|}{t})
        = - \frac{\left|\mathbf{q}\right|}{t^2} \nabla{w(\frac{\left|\mathbf{q}\right|}{t})} \f$

        where \f$ \left|\mathbf{q}\right| \f$ represents the norm of the
        query coordinates expressed in centered basis.
    */
    PONCA_MULTIARCH inline Scalar scaledw(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

    /*!
        \brief Second order derivative in scale  \f$t\f$

        \f$ \frac{\delta^2 \frac{\left|\mathbf{q}\right|}{t}}{\delta t^2}
        \nabla w(\frac{\left|\mathbf{q}\right|}{t}) +
        \left(\frac{\delta \frac{\left|\mathbf{q}\right|}{t}}{\delta t}\right)^2
        \nabla^2 w(\frac{\left|\mathbf{q}\right|}{t}) =
        \frac{2\left|\mathbf{q}\right|}{t^3} \nabla{w(\frac{\left|\mathbf{q}\right|}{t})} +
        \frac{\left|\mathbf{q}\right|^2}{t^4} \nabla^2{w(\frac{\left|\mathbf{q}\right|}{t})}
        \f$

        where \f$ \left|\mathbf{q}\right| \f$ represents the norm of the
        query coordinates expressed in centered basis.
    */
    PONCA_MULTIARCH inline Scalar scaled2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

    /*!
        \brief Cross derivative in scale \f$t\f$ and in space (for each spatial dimension \f$\mathsf{x})\f$

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
    */
    PONCA_MULTIARCH inline VectorType scaleSpaced2w(const VectorType& _q,
        const DataPoint&  /*attributes*/) const;

    /*! \brief Access to the evaluation scale set during the initialization */
    PONCA_MULTIARCH inline Scalar evalScale() const { return m_t; }

    /*! \brief Access to the evaluation position set during the initialization */
    PONCA_MULTIARCH inline const VectorType & evalPos() const { return m_p; }

protected:
    Scalar       m_t;  /*!< \brief Evaluation scale */
    WeightKernel m_wk; /*!< \brief 1D function applied to weight queries */
    VectorType   m_p;  /*!< \brief basis center */

};// class DistWeightFunc

#include "weightFunc.hpp"

}// namespace Ponca
