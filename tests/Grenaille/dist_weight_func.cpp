/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/dist_weight_func.cpp
    \brief Test distance weight function derivatives
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename WeightKernel>
void testFunction()
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;

    DataPoint dummy;

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar t       = Eigen::internal::random<Scalar>(0.10, 10.0);
    Scalar h       = Scalar(1e-5);
    Scalar h2      = Scalar(.5)*h;

    DistWeightFunc<DataPoint, WeightKernel> wfunc(t);
    DistWeightFunc<DataPoint, WeightKernel> wfuncR(t+h);
    DistWeightFunc<DataPoint, WeightKernel> wfuncL(t-h);
    DistWeightFunc<DataPoint, WeightKernel> wfuncR2(t+h2);
    DistWeightFunc<DataPoint, WeightKernel> wfuncL2(t-h2);

    // rejection sampling inside a sphere of radius t
    VectorType x = t*VectorType::Random();
    while ( x.norm() < 0.001*t || 0.999*t < x.norm())
        x = t*VectorType::Random();

    // actual values of derivatives
    VectorType dx_w   = wfunc.spacedw(x, dummy);
    MatrixType d2x_w  = wfunc.spaced2w(x, dummy);
    Scalar     dt_w   = wfunc.scaledw(x, dummy);
    Scalar     d2t_w  = wfunc.scaled2w(x, dummy);
    VectorType d2tx_w = wfunc.scaleSpaced2w(x, dummy);

    // each derivative is approximated by a centerd finit difference

    // 1st order spatial derivative approximation
    VectorType dx_w_;
    for(int i=0; i<DataPoint::Dim; ++i)
    {
        VectorType ei = VectorType::Zero();
        ei[i] = Scalar(1.);
        dx_w_[i] = (wfunc.w(x+h*ei, dummy) - wfunc.w(x-h*ei, dummy))/(Scalar(2.)*h);
    }

    // 2nd order spatial derivative approximation
    MatrixType d2x_w_;
    for(int i=0; i<DataPoint::Dim; ++i)
    {
        VectorType ei = VectorType::Zero();
        ei[i] = Scalar(1.);
        for(int j=0; j<DataPoint::Dim; ++j)
        {
            VectorType ej = VectorType::Zero();
            ej[j] = Scalar(1.);

            d2x_w_(i,j) = Scalar(1.)/(h*h)*(
                            + wfunc.w(x + h2*ei + h2*ej, dummy)
                            - wfunc.w(x + h2*ei - h2*ej, dummy)
                            - wfunc.w(x - h2*ei + h2*ej, dummy)
                            + wfunc.w(x - h2*ei - h2*ej, dummy)
                        );
        }
    }

    // 1st order scale derivative approximation
    Scalar dt_w_ = (wfuncR.w(x, dummy) - wfuncL.w(x, dummy))/(Scalar(2.)*h);

    // 2nd order scale derivative approximation
    Scalar d2t_w_ = (wfuncR.w(x, dummy)
                     -Scalar(2.)*wfunc.w(x, dummy)
                     + wfuncL.w(x, dummy))/(h*h);

    // cross derivatives approximation
    VectorType d2tx_w_;
    for(int i=0; i<DataPoint::Dim; ++i)
    {
        VectorType ei = VectorType::Zero();
        ei[i] = Scalar(1.);
        d2tx_w_[i] = Scalar(1.)/(h*h)*(
          + wfuncR2.w(x + h2*ei, dummy)
          - wfuncR2.w(x - h2*ei, dummy)
          - wfuncL2.w(x + h2*ei, dummy)
          + wfuncL2.w(x - h2*ei, dummy)
        );
    }

    for(int i=0; i<DataPoint::Dim; ++i)
    {
        VERIFY( std::abs(dx_w[i]-dx_w_[i]) < epsilon);
        VERIFY( std::abs(d2tx_w[i]-d2tx_w_[i]) < epsilon);
        for(int j=0; j<DataPoint::Dim; ++j)
        {
            VERIFY( std::abs(d2x_w(i,j)-d2x_w_(i,j)) < epsilon);
        }
    }
    VERIFY( std::abs(dt_w-dt_w_) < epsilon );
    VERIFY( std::abs(d2t_w-d2t_w_) < epsilon );
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> DataPoint;

    typedef SmoothWeightKernel<Scalar> SmoothKernel;
    typedef ConstantWeightKernel<Scalar> ConstantKernel;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<DataPoint, SmoothKernel>() ));
        CALL_SUBTEST(( testFunction<DataPoint, ConstantKernel>() ));
    }
    cout << "ok" << endl;
}


int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Verify distance weight function derivatives" << endl;

    callSubTests<long double, 1>();
    callSubTests<long double, 2>();
    callSubTests<long double, 3>();
    callSubTests<long double, 4>();
    callSubTests<long double, 5>();
}
