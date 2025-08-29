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
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <chrono>
#include <math.h>
#include <iomanip>

using namespace std;
using namespace Ponca;

/*
 * Test the OrientedSphereFit on a paraboloid using a given neighborhood filter
 */
template<typename DataPoint, typename NF>
void testFunction(typename DataPoint::Scalar lowPrecisionEpsilon = typename DataPoint::Scalar(0.001)) // Lesser precision for the paraboloid test
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef Basket<DataPoint, NF, OrientedSphereFit> Fit;

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
    VectorType center = 1000 * VectorType::Random();

    Scalar zmax = std::abs((coeff[0] + coeff[1]) * width*width);
    Scalar analysisScale = std::sqrt(zmax*zmax + width*width);

    Fit fit;
    fit.setNeighborFilter(NF(center, analysisScale));
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

    if(fit.isStable())
    {
        std::vector<VectorType> samples (nbPoints);
        for (int i = 0; i < nbPoints; ++i)
        {
            VectorType p = center + analysisScale * VectorType::Random();
            samples[i] = p;

            // check that the projected point is on the surface
            VERIFY( std::abs(fit.potential(fit.project(p))) < lowPrecisionEpsilon );
        }

        // check if direct projection gives same or better result than descent projection.
        for( const auto& p: samples )
        {
            VectorType res1 = fit.project( p );
            VectorType res2 = fit.projectDescent( p, 1000 ); // force high number of iterations
            VERIFY( res1.isApprox( res2, lowPrecisionEpsilon ));
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

/*
 * Compare two neighborhood filter using the OrientedSphereFit on the paraboloid test
 */
template<typename DataPoint, typename NF, typename NFG>
void compareNF(typename DataPoint::Scalar lowPrecisionEpsilon = typename DataPoint::Scalar(0.1)) // Lesser precision for the paraboloid test
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef Basket<DataPoint, NF, OrientedSphereFit> Fit;
    typedef Basket<DataPoint, NFG, OrientedSphereFit> FitG;

    // generate samples
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
    VectorType center = 1000 * VectorType::Random();

    Scalar zmax = std::abs((coeff[0] + coeff[1]) * width*width);
    Scalar analysisScale = std::sqrt(zmax*zmax + width*width);

    Fit fit;
    FitG fitG;
    fit.setNeighborFilter(NF(center, analysisScale));
    fitG.setNeighborFilter(NFG(center, analysisScale));
    fit.init();
    fitG.init();

    for(int i = 0; i < nbPointsFit; ++i)
    {
        DataPoint p = getPointOnParaboloid<DataPoint>(coeff.x(),
                                                      coeff.y(),
                                                      width,
                                                      false);           // noise
        p.pos() += center;

        fit.addNeighbor(p);
        fitG.addNeighbor(DataPoint(p.pos(), p.normal()));
    }
    std::cout << "------------------------------------------------ IN fit.finalize() ------------------------------------------------ " << std::endl;
    fit.finalize();
    std::cout << "------------------------------------------------ IN fitG.finalize() ----------------------------------------------- " << std::endl;
    fitG.finalize();
    std::cout << "------------------------------------------------------------------------------------------------------------------- " << std::endl;

    if(fit.isStable() && fitG.isStable())
    {
        std::vector<VectorType> samples (nbPoints);
        for (int i = 0; i < nbPoints; ++i)
        {
            VectorType p = center + analysisScale * VectorType::Random();
            samples[i] = p;

            if (std::abs(fitG.potential(fitG.project(p))) >= lowPrecisionEpsilon)
            {
                std::cout << "------------------------------------------------ IN testFunction() ------------------------------------------------ " << std::endl;
            	std::cout << std::left << std::setw(26);
                std::cout << "fit.project(p) : " << fit.project(p).transpose();
            	std::cout << std::right << std::setw(25);
                std::cout << "| fitG.project(p) : " << fitG.project(p).transpose();
                std::cout << "\n";

            	std::cout << std::left << std::setw(25);
                std::cout << "fit.potential(fit.project(p)) : " << fit.potential(fit.project(p));
            	std::cout << std::right << std::setw(60);
                std::cout << "| fitG.potential(fitG.project(p)) : " << fitG.potential(fitG.project(p));
                std::cout << "\n" << endl;

            }

            // check that the projected point is on the surface
            VERIFY( std::abs(fit.potential(fit.project(p))) < lowPrecisionEpsilon );
            VERIFY( std::abs(fitG.potential(fitG.project(p))) < lowPrecisionEpsilon );
        }

        // check if direct projection gives same or better result than descent projection.
        for( const auto& p: samples )
        {
            VectorType res1 = fit.project( p );
            VectorType res2 = fit.projectDescent( p, 1000 ); // force high number of iterations
            VERIFY( res1.isApprox( res2, lowPrecisionEpsilon ));
            VectorType res1G = fitG.project( p );
            VectorType res2G = fitG.projectDescent( p, 1000 ); // force high number of iterations
            VERIFY( res1G.isApprox( res2G, lowPrecisionEpsilon ));
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;
    typedef DistWeightFuncGlobal<Point, ConstantWeightKernel<Scalar> > WeightConstantFuncGlobal;
    typedef DistWeightFuncGlobal<Point, SmoothWeightKernel<Scalar> > WeightSmoothFuncGlobal;

    cout << "Testing with parabola..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( compareNF<Point, WeightSmoothFunc, WeightSmoothFuncGlobal>(0.1) ));
        CALL_SUBTEST(( compareNF<Point, WeightConstantFunc, WeightConstantFuncGlobal>(0.1) ));
        // CALL_SUBTEST(( testFunction<Point, WeightSmoothFuncGlobal>() ));
        // CALL_SUBTEST(( testFunction<Point, WeightConstantFuncGlobal>() ));
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
