/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/smooth_weight_kernel.cpp
    \brief Test smooth weight kernel derivatives
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

using namespace std;
using namespace Grenaille;

template<class Kernel>
void testFunction()
{
    typedef typename Kernel::Scalar Scalar;

    Scalar step = Scalar(0.05);
    int n = Scalar(1.)/step;

    Kernel k;
    Scalar h = Scalar(1e-6);
    Scalar epsilon = testEpsilon<Scalar>();

    // compare to finite difference approximations
    for(int i=1; i<=n; ++i)
    {
        Scalar a = i*step;

        Scalar f    = k.f(a);
        Scalar fr   = k.f(a+h);
        Scalar fl   = k.f(a-h);

        Scalar df   = k.df(a);
        Scalar ddf  = k.ddf(a);

        Scalar df_  = (fr - fl)/(Scalar(2.)*h);
        Scalar ddf_ = (fr - Scalar(2.)*f + fl)/(h*h);

        Scalar diff1 = std::abs(df-df_);
        Scalar diff2 = std::abs(ddf-ddf_);

        VERIFY(diff1 < epsilon);
        VERIFY(diff2 < epsilon);
    }
}

template<typename Scalar>
void callSubTests()
{
    typedef SmoothWeightKernel<Scalar> Kernel;
    CALL_SUBTEST(( testFunction<Kernel>() ));
    cout << "ok" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Verify smooth weight kernel derivatives" << endl;

    callSubTests<long double>();
}
