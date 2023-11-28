/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/capability.cpp
    \brief Test optional capabilities
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/capabilities.h>
#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/plane.h>
#include <Ponca/src/Fitting/algebraicSphere.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

using namespace Ponca;

constexpr int Dim = 3;
using DataPoint = PointPosition<float, Dim>;
using Kernel = ConstantWeightKernel<float>;
using WeightFunc = DistWeightFunc<DataPoint, Kernel>;

// compile-time and run-time check if a given value is true or false
// this template function is required because "Outside a template, a discarded 
// statement is fully checked." (cppreference.com/w/cpp/language/if) 
template<bool b>
void CHECK() {
    static_assert(b);
    VERIFY(b);
}

// extension that test if PLANE is provided
template<class P, class W, typename T>
class InitTestExtensionPlane : public T
{
public:
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    
    inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        Base::init(_basisCenter);
        if constexpr( IS_PROVIDED(PLANE) ) {
            CHECK<true>();
            Base::setPlane(_basisCenter, VectorType::UnitZ());
        } else {
            // must not compile if PLANE is provided
            CHECK<false>();
        }

        if constexpr( IS_PROVIDED(DOES_NOT_EXIST) ) {
            // must not compile because DOES_NOT_EXIST does not exist
            CHECK<false>();
        }
    }
};

// extension that test if ALGEBRAIC_SPHERE is provided
template<class P, class W, typename T>
class InitTestExtensionSphere : public T
{
public:
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    
    inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        Base::init(_basisCenter);
        if constexpr( IS_PROVIDED(ALGEBRAIC_SPHERE) ) {
            CHECK<true>();
            Base::m_uc = 0;
            Base::m_ul = VectorType::UnitZ();
            Base::m_uq = 0;
        } else {
            // must not compile if ALGEBRAIC_SPHERE is provided
            CHECK<false>();
        }

        if constexpr( IS_PROVIDED(DOES_NOT_EXIST) ) {
            // must not compile because DOES_NOT_EXIST does not exist
            CHECK<false>();
        }
    }
};

using PlaneFit  = Basket<DataPoint, WeightFunc, Plane, InitTestExtensionPlane>;
using SphereFit = Basket<DataPoint, WeightFunc, AlgebraicSphere, InitTestExtensionSphere>;

using ErrorFit1 = Basket<DataPoint, WeightFunc, Plane, InitTestExtensionSphere>;
using ErrorFit2 = Basket<DataPoint, WeightFunc, AlgebraicSphere, InitTestExtensionPlane>;

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    PlaneFit plane;
    plane.init();
    
    SphereFit sphere;
    sphere.init();

    ErrorFit1 error1;
    // error1.init();   // error: static assertion failed (because ALGEBRAIC_SPHERE is not provided)

    ErrorFit2 error2;
    // error2.init();   // error: static assertion failed (because PLANE is not provided)

    return EXIT_SUCCESS;
}
