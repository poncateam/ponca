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
    \brief Plane fitting procedure on oriented point sets using robust implicit mean local surfaces

    Method published in \cite oztireli:2009:feature

    A plane is fitted by iteratively minimizing the following equation

    \f$ f^k(x) = argmin_{S_0} \sum ( S_0 + ( x_i - x )^T n_i )^2 &Phi; _i(x) w(r_i^{k-1} ) 
    w_n( &Delta; n_i^k )\f$

    with the residuals \f$ r_i^{k-1} = f^{k-1}(x) - ( x - x_i )^T n_i \f$

    \f$w\f$ and \f$w_n\f$ are Gaussian functions

    &Phi; is an arbitrary weighting function

    \inherit Concept::FittingProcedureConcept

    \see Plane

 */
template<class DataPoint, class _WFunctor, typename T>
class RimlsPlaneFitImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum { Check = Base::PROVIDES_PLANE };

protected:
    /*!< \brief A parameter which controls the sharpness of the result, smaller value lead to sharper result */
    Scalar m_sigmaN;

    /*!< \brief convergence threshold ( default = \f$ 10^{-3} \f$ ) */
    Scalar m_minConvergence;

    /*!< \brief maximum number of iteration  ( default = \f$ 10 \f$ ) */
    unsigned int m_maxIteration;

    /*!< \brief neighbors are accumulated in a vector */
    std::vector<DataPoint> m_neighbors;

public:
    PONCA_MULTIARCH inline RimlsPlaneFitImpl() : Base() {}

    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    PONCA_MULTIARCH inline const Scalar sigmaN() const { return m_sigmaN; }
    PONCA_MULTIARCH inline Scalar& sigmaN() { return m_sigmaN; }

    PONCA_MULTIARCH inline const Scalar minConvergence() const { return m_minConvergence; }
    PONCA_MULTIARCH inline Scalar& minConvergence() { return m_minConvergence; }

    PONCA_MULTIARCH inline const unsigned int maxIteration() const { return m_maxIteration; }
    PONCA_MULTIARCH unsigned int& maxIteration() { return m_maxIteration; }

private:
    PONCA_MULTIARCH inline Scalar gaussian(Scalar x, Scalar sigma)
    {
        PONCA_MULTIARCH_STD_MATH(exp);
        return exp(-(x) / (sigma * sigma));
    }

}; // class RimlsPlaneFitImpl

/// @brief Helper alias for rimls plane fitting on points using RimlsPlaneFittingImpl
//! [RimlsPlaneFit Definition]
    template< class DataPoint, class _WFunctor, typename T >
    using RimlsPlaneFit = 
        RimlsPlaneFitImpl< DataPoint, _WFunctor, 
        Plane< DataPoint, _WFunctor, T > >;
//! [RimlsPlaneFit Definition]

#include "rimlsPlaneFit.hpp"
} //namespace Ponca

