/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/weight_kernel.cpp
    \brief Test weight kernel derivatives
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <unsupported/Eigen/AutoDiff>

#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

using namespace std;
using namespace Ponca;

template<class Kernel>
void testFunctionAutoDiff()
{
    typedef typename Kernel::Scalar ScalarDiff;
    typedef typename ScalarDiff::Scalar Scalar;

    Scalar step = Scalar(0.05);
    int n = Scalar(1.)/step;

    Kernel k;
    Scalar epsilon = testEpsilon<Scalar>();

    // compare to automatic differentiation
    for(int i=1; i<=n; ++i)
    {
        ScalarDiff a(i*step, 1, 0);

        ScalarDiff f   = k.f(a);
        ScalarDiff df  = k.df(a);
        ScalarDiff ddf = k.ddf(a);

        Scalar diff1 = std::abs( f.derivatives()[0] - df.value());
        Scalar diff2 = std::abs(df.derivatives()[0] - ddf.value());

        VERIFY(diff1 < epsilon);
        VERIFY(diff2 < epsilon);
    }
}

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

template<typename Scalar, template <typename > class KernelT>
void callSubTests()
{
    typedef Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,1,1>> ScalarDiff;

    typedef KernelT<Scalar> Kernel;
    typedef KernelT<ScalarDiff> KernelAutoDiff;

    CALL_SUBTEST(( testFunction<Kernel>() ));
    CALL_SUBTEST(( testFunctionAutoDiff<KernelAutoDiff>() ));

    cout << "ok" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Verify smooth weight kernel derivatives" << endl;

    callSubTests<long double, SmoothWeightKernel>();

    cout << "Verify Wendland weight kernel derivatives" << endl;

    callSubTests<long double, WendlandWeightKernel>();

    cout << "Verify singular weight kernel derivatives" << endl;

    callSubTests<long double, SingularWeightKernel>();
}
