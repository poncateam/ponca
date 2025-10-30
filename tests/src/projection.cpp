/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble, Nicolas Mellado
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/sphereFit.h>
#include <Ponca/src/Fitting/unorientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <chrono>
#include <math.h>

using namespace std;
using namespace Ponca;

/*
 * Test the OrientedSphereFit on a paraboloid using a given neighborhood filter
 */
template<typename Fit>
void testFunction(typename Fit::Scalar lowPrecisionEpsilon = typename Fit::Scalar(0.001)) // Lesser precision for the paraboloid test
{
    // Define related structure
    typedef typename Fit::Scalar Scalar;
    typedef typename Fit::VectorType VectorType;
    typedef typename Fit::DataPoint DataPoint;

    //generate samples
    int nbPoints = Eigen::internal::random<int>(100, 1000);
    int nbPointsFit = 50;

    // equal probability of having a plane or a random quadric
    VectorType coeff = 5 * VectorType::Random();
    if(Eigen::internal::random<Scalar>(0., 1.) < Scalar(0.5))
    {
        coeff = VectorType::Zero();
    }
    // avoid saddles, which can cause issues with the sphere fitting
    coeff.y() = std::copysign(coeff.y(), coeff.x());

    Scalar width = Eigen::internal::random<Scalar>(1., 10.);
    // maximum offset is <5 unit, and always smaller than 2*width.
    // It is plenty enough to test local/global basis robustness, without introducing rounding errors
    Scalar offset = Eigen::internal::random<Scalar>(1., std::min(Scalar(5.), Scalar(.5)*width));
    VectorType center = offset*VectorType::Random();

    Scalar zmax = std::abs((coeff[0] + coeff[1]) * width*width);
    Scalar analysisScale = std::sqrt(zmax*zmax + width*width + offset);

    Fit fit;
    fit.setNeighborFilter({center, analysisScale});
    fit.init();

    for(int i = 0; i < nbPointsFit; ++i)
    {
        DataPoint p = getPointOnParaboloid<DataPoint>(coeff.x(),
                                                      coeff.y(),
                                                      width,
                                                      false);           // noise
        p.pos() += center;

        fit.addNeighbor(p);
    }
    fit.finalize();

    if(fit.isStable())
    {
        for (int i = 0; i < nbPoints; ++i)
        {
            const VectorType p = center + analysisScale * VectorType::Random();

            // check that the projected point is on the surface
            VectorType projD = fit.projectDescent( p, 1000 );
            VERIFY( std::abs(fit.potential(projD)) < lowPrecisionEpsilon );

            VectorType proj = fit.project( p );
            Scalar p1 = std::abs(fit.potential(proj));
            Scalar p2 = std::abs(fit.potential(p));
            // check the direct projection did not move the point away from the surface (can be stationary if already
            // on the surface)
            VERIFY( p1 <= p2 );
        }

        // Disable this test: not true with apple-clang 12.
#ifdef COMPARE_PROJECTION_TIMINGS
        auto start1 = std::chrono::system_clock::now();
        for( const auto& p: samples )
          fit.project( p );
        auto end1 = std::chrono::system_clock::now();

        auto start2 = std::chrono::system_clock::now();
        for( const auto& p: samples )
          fit.projectDescent( p );
        auto end2 = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds1 = end1-start1;
        std::chrono::duration<double> elapsed_seconds2 = end2-start2;
        std::cout << "Default: " << elapsed_seconds1.count() << " Descent: " << elapsed_seconds2.count() << "s\n";
        VERIFY( elapsed_seconds1 <= elapsed_seconds2 );
#endif
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    using Point = PointPositionNormal<Scalar, Dim> ;

    using WeightSmoothFunc        = DistWeightFunc<Point, SmoothWeightKernel<Scalar>>;
    using WeightConstantFuncLocal = Ponca::DistWeightFunc<Point, Ponca::ConstantWeightKernel<Scalar>>;
    using NoWeightFuncGlobal      = Ponca::NoWeightFuncGlobal<Point> ;
    using NoWeightFunc            = Ponca::NoWeightFunc<Point> ;

#define MAKE_FIT_TYPE(Fit,Weight) Basket<Point, Weight, Fit>

#define TEST_FIT(Fit) \
        CALL_SUBTEST(( testFunction<MAKE_FIT_TYPE(Fit,WeightSmoothFunc)>() )); \
        CALL_SUBTEST(( testFunction<MAKE_FIT_TYPE(Fit,WeightConstantFuncLocal)>() )); \
        CALL_SUBTEST(( testFunction<MAKE_FIT_TYPE(Fit,NoWeightFunc)>() ));

    cout << "Testing with " << typeid(Scalar).name() << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        TEST_FIT(OrientedSphereFit) // AlgebraicSphere requires local basis
        TEST_FIT(CovariancePlaneFit)
        CALL_SUBTEST(( testFunction<MAKE_FIT_TYPE(CovariancePlaneFit,NoWeightFuncGlobal)>() ));
        TEST_FIT(UnorientedSphereFit)
        TEST_FIT(SphereFit)

    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test projection for different baskets..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
