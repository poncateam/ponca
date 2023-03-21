/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

/*!
    \file weightKernel.h Define 1D weight kernel functors
*/


namespace Ponca
{
/*!
    \brief Concept::WeightKernelConcept returning a constant value

    \inherit Concept::WeightKernelConcept
*/
template <typename _Scalar>
class ConstantWeightKernel
{
public:
    /*! \brief Scalar type defined outside the class */
    typedef _Scalar Scalar;

    // Init
    //! \brief Default constructor that could be used to set the returned value
    PONCA_MULTIARCH inline ConstantWeightKernel(const Scalar& _value = Scalar(1.)) : m_y(_value){}
    //! \brief Set the returned value
    PONCA_MULTIARCH inline void setValue(const Scalar& _value){ m_y = _value; }

    // Functor
    //! \brief Return the constant value
    PONCA_MULTIARCH inline Scalar f  (const Scalar&) const { return m_y; }
    //! \brief Return \f$ 0 \f$
    PONCA_MULTIARCH inline Scalar df (const Scalar&) const { return Scalar(0.); }
    //! \brief Return \f$ 0 \f$
    PONCA_MULTIARCH inline Scalar ddf(const Scalar&) const { return Scalar(0.); }

private:
    Scalar m_y; /*!< \brief Constant value returned by the kernel */
};// class ConstantWeightKernel


/*!
    \brief Smooth WeightKernel defined in \f$\left[0 : 1\right]\f$
    \todo Add a degree value as template parameter (in this class or another one), with specialized functions for 2

    \inherit Concept::WeightKernelConcept
*/
template <typename _Scalar>
class SmoothWeightKernel
{
public:
    /*! \brief Scalar type defined outside the class*/
    typedef _Scalar Scalar;

    // Functor
    /*! \brief Defines the smooth weighting function \f$ w(x) = (x^2-1)^2 \f$ */
    PONCA_MULTIARCH inline Scalar f  (const Scalar& _x) const { Scalar v = _x*_x - Scalar(1.); return v*v; }
    /*! \brief Defines the smooth first order weighting function \f$ \nabla w(x) = 4x(x^2-1) \f$ */
    PONCA_MULTIARCH inline Scalar df (const Scalar& _x) const { return Scalar(4.)*_x*(_x*_x-Scalar(1.)); }
    /*! \brief Defines the smooth second order weighting function \f$ \nabla^2 w(x) = 12x^2-4 \f$ */
    PONCA_MULTIARCH inline Scalar ddf(const Scalar& _x) const { return Scalar(12.)*_x*_x - Scalar(4.); }
};//class SmoothWeightKernel


/*!
    \brief Wendland WeightKernel defined in \f$\left[0 : 1\right]\f$

    \inherit Concept::WeightKernelConcept

    Weight function is an implementation of equation 2 in \cite Alexa:2009:Hermite
*/

template <typename _Scalar>
class WendlandWeightKernel
{
public:
    /*! \brief Scalar type defined outside the class*/
    typedef _Scalar Scalar;

    // Functor
    /*! \brief Defines the Wendland weighting function \f$ w(x) = (1-x)^4(4x+1) \f$ */
    PONCA_MULTIARCH inline Scalar f  (const Scalar& _x) const {
        const Scalar v = Scalar(1.) - _x;
        return v * v * v * v * ((Scalar(4.) * _x) + Scalar(1.));
    }
    /*! \brief Defines the Wendland first order weighting function \f$ \nabla w(x) = 20x * (x−1)^3 \f$ */
    PONCA_MULTIARCH inline Scalar df (const Scalar& _x) const {
        const Scalar v = _x - Scalar(1.);
        return Scalar(20.) * _x * v * v * v;
    }
    /*! \brief Defines the Wendland second order weighting function \f$ \nabla^2 w(x) = (x−1)^2 * (80x−20) \f$ */
    PONCA_MULTIARCH inline Scalar ddf(const Scalar& _x) const {
        const Scalar v = _x - Scalar(1.);
        return v * v * (Scalar(80.) * _x - Scalar(20));
    }
};//class WendlandWeightKernel


/*!
    \brief Singular WeightKernel defined in \f$\left[0 : 1\right]\f$

    \inherit Concept::WeightKernelConcept

    Weight function is an implementation of an unnumbered equation but defined in the Appendices in \cite Alexa:2009:Hermite
*/

template <typename _Scalar>
class SingularWeightKernel
{
public:
    /*! \brief Scalar type defined outside the class*/
    typedef _Scalar Scalar;

    // Functor
    /*! \brief Defines the Singular weighting function \f$ w(x) = 1 / (x^2) \f$ */
    PONCA_MULTIARCH inline Scalar f  (const Scalar& _x) const {
        return Scalar(1.) / (_x * _x);
    }
    /*! \brief Defines the Singular first order weighting function \f$ \nabla w(x) = -2 / (x^3) \f$ */
    PONCA_MULTIARCH inline Scalar df (const Scalar& _x) const {
        return Scalar(-2.) / (_x * _x * _x);
    }
    /*! \brief Defines the Singular second order weighting function \f$ \nabla^2 w(x) = 6 / (x^4) \f$ */
    PONCA_MULTIARCH inline Scalar ddf(const Scalar& _x) const {
        return Scalar(6.) / (_x * _x * _x * _x);
    }
};//class SingularWeightKernel


/*!
    \brief Compact Exponential WeightKernel defined in \f$\left[0 : 1\right]\f$

    Continuity is \f$C^\infty\f$ for \f$x=1\f$.

    \warning \f$\nabla^2 w(x)\f$ is not defined for all values in \f$\left[0 : 1\right]\f$. Do not use for second-order differentiation.

    \inherit Concept::WeightKernelConcept
*/
template <typename _Scalar>
class CompactExpWeightKernel
{
public:
    /*! \brief Scalar type defined outside the class*/
    typedef _Scalar Scalar;

    // Functor
    /*! \brief Defines the smooth weighting function \f$ w(x) = e^{-\frac{x^2}{1 - x^2}} \f$
     *  \see https://www.wolframalpha.com/input?i=e%5E%28-x%5E2%2F%281+-+x%5E2%29%29&assumption=%22ClashPrefs%22+-%3E+%7B%22Math%22%7D
     */
    PONCA_MULTIARCH inline Scalar f  (const Scalar& _x) const { Scalar v = _x*_x; return exp(-v/(Scalar(1)-v)); }
    /*! \brief Defines the smooth first order weighting function \f$ \nabla w(x) = -\frac{2 x e^{\frac{x^2}{x^2 - 1}}}{(1 - x^2)^2} \f$
     * \see https://www.wolframalpha.com/input?i2d=true&i=+-Divide%5B%5C%2840%292+Power%5Be%2C%5C%2840%29Power%5Bx%2CDivide%5B2%2C%5C%2840%29Power%5Bx%2C2%5D+-+1%5C%2841%29%5D%5D%5C%2841%29%5D+x%5C%2841%29%2CPower%5B%5C%2840%291+-+Power%5Bx%2C2%5D%5C%2841%29%2C2%5D%5D
     */
    PONCA_MULTIARCH inline Scalar df (const Scalar& _x) const { Scalar v = _x*_x; Scalar mv = v-Scalar(1); return -(Scalar(2) * exp(v/mv) * _x)/(mv*mv); }
    /*! \brief Defines the smooth second order weighting function
     * \f$ \nabla^2 w(x) = \frac{2 e^\frac{x^2}{x^2 - 1} \left(4 x^{\frac{2}{x^2 - 1} + 2} log(x) - (x^2 - 1) \left(-3 x^2 + 2 x^\frac{2}{x^2 - 1} - 1\right)\right)}{(x^2 - 1)^4} \f$
     * \see \see https://www.wolframalpha.com/input?i2d=true&i=+-Divide%5B%5C%2840%292+Power%5Be%2C%5C%2840%29Power%5Bx%2CDivide%5B2%2C%5C%2840%29Power%5Bx%2C2%5D+-+1%5C%2841%29%5D%5D%5C%2841%29%5D+x%5C%2841%29%2CPower%5B%5C%2840%291+-+Power%5Bx%2C2%5D%5C%2841%29%2C2%5D%5D
     */
    PONCA_MULTIARCH inline Scalar ddf(const Scalar& _x) const {
        Scalar v = _x*_x;            // x^2
        Scalar mv = v-Scalar(1);     // x^2-1
        Scalar mvpow = Scalar(2)/mv; // 2/(x^2-1)
        Scalar cmvpow = pow(_x,mvpow); // x^{2/(x^2-1)}
        return Scalar(2) * exp(cmvpow) * ((Scalar(4) * pow(_x,mvpow + Scalar(2)) * log(_x) - mv * (-Scalar(3) * v + Scalar(2) * cmvpow - Scalar(1))))/(mv*mv*mv*mv);
    }

    //! \brief #df is defined and valid on the definition interval
    static constexpr bool isDValid = true;
    //! \brief #ddf is not defined and valid on the definition interval
    static constexpr bool isDDValid = false;
};//class CompactExpWeightKernel


}// namespace Ponca
