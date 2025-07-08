/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/Grenaille/okabe_primitive.cpp
    \brief Test validity Algebraic Sphere Primitive
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>


using namespace std;
using namespace Ponca;

template<typename Point, typename Fit, typename WeightFunc>
void testFunction()
{
        // Define related structure
    typedef typename Point::Scalar Scalar;
    typedef typename Point::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);

    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,1000);
    VectorType center = VectorType::Random() * centerScale;

    Scalar epsilon = testEpsilon<Scalar>();

    vector<Point> vecs ( nbPoints );

    for(unsigned int i = 0; i < vecs.size(); ++i)
    {
        vecs[i] = getPointOnSphere<Point>(radius, center, true, false, false);
    }


#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int k=0; k<int(vecs.size()); ++k)
    {
        const auto &fitInitPos = vecs[k].pos();

        Fit fit;
        fit.setWeightFunc(WeightFunc(fitInitPos, analysisScale));
        fit.compute(vecs);

        if(fit.isStable())
        {
            Fit f2 = fit;
            // do two iterations to test successive calls
            VectorType candidate = VectorType::Random();
            for (unsigned int j = 0; j != 2; ++j){
                // move basis center to a portion of the radius of the sphere
                f2.changeBasis(fitInitPos + Scalar(0.1)*fit.radius()*VectorType::Random());

                VERIFY(  fit.radius() - f2.radius() < epsilon );
                VERIFY( (fit.center() - f2.center() ).norm() < epsilon );
                VERIFY( (fit.project(candidate) - f2.project(candidate)).norm() - Scalar(1.) < epsilon );
            }
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    using Point = PointPositionNormal<Scalar, Dim>;

    // We test only primitive functions and not the fitting procedure
    using WeightFunc = DistWeightFunc<Point, SmoothWeightKernel<Scalar> >;
    using Sphere     = Basket<Point, WeightFunc, OrientedSphereFit>;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, Sphere, WeightFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Algebraic Sphere Primitive functions in 3 dimensions..." << endl;
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();

    cout << "Test Algebraic Sphere Primitive functions in 5 dimensions..." << endl;
    callSubTests<float, 5>();
    callSubTests<double, 5>();
    callSubTests<long double, 5>();
    cout << "Ok..." << endl;
}
