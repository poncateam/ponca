/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/fit_unoriented.cpp
 \brief Test coherance of unoriented sphere fitting
 */


#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/gls.h>
#include <Ponca/src/Fitting/unorientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(Scalar(0.1), Scalar(10.));

    Scalar analysisScale = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10);
    VectorType center = VectorType::Random() * centerScale;

    Scalar epsilon = testEpsilon<Scalar>();
    // epsilon is relative to the radius size
//    Scalar radiusEpsilon = epsilon * radius;


    vector<DataPoint> vectorPoints(nbPoints);
    vector<DataPoint> vectorReversedNormals100(nbPoints);
    vector<DataPoint> vectorReversedNormalsRandom(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        //reverse random normals
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise);
    }

    reverseNormals<DataPoint>(vectorReversedNormals100, vectorPoints, false);
    reverseNormals<DataPoint>(vectorReversedNormalsRandom, vectorPoints);

    // Test sphere descriptors coherance for each points
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit, fitReverse100, fitReverseRandom;

        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints);

        fitReverse100.setWeightFunc(WeightFunc(analysisScale));
        fitReverse100.init(vectorReversedNormals100[i].pos());
        fitReverse100.compute(vectorPoints);

        fitReverseRandom.setWeightFunc(WeightFunc(analysisScale));
        fitReverseRandom.init(vectorReversedNormalsRandom[i].pos());
        fitReverseRandom.compute(vectorPoints);

        if(fit.isStable() && fitReverse100.isStable() && fitReverseRandom.isStable())
        {
            Scalar kappa1 = std::abs(fit.kappa());
            Scalar kappa2 = std::abs(fitReverse100.kappa());
            Scalar kappa3 = std::abs(fitReverseRandom.kappa());

            Scalar tau1 = std::abs(fit.tau());
            Scalar tau2 = std::abs(fitReverse100.tau());
            Scalar tau3 = std::abs(fitReverseRandom.tau());

            VectorType eta1 = fit.eta().normalized().array().abs();
            VectorType eta2 = fitReverse100.eta().normalized().array().abs();
            VectorType eta3 = fitReverseRandom.eta().normalized().array().abs();

            // Check kappa coherance
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(kappa1 - kappa2), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(kappa1 - kappa3), Scalar(1.), epsilon) );

            // Check tau coherance
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(tau1 - tau2), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(tau1 - tau3), Scalar(1.), epsilon) );

            // Check eta coherance
            VERIFY( Eigen::internal::isMuchSmallerThan((eta1 - eta2).norm(), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan((eta1 - eta3).norm(), Scalar(1.), epsilon) );
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, UnorientedSphereFit, GLSParam> FitSmoothUnoriented;
    typedef Basket<Point, WeightConstantFunc, UnorientedSphereFit, GLSParam> FitConstantUnoriented;

    cout << "Testing with perfect sphere (unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect sphere
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test unoriented sphere fit coherance..." << endl;

    callSubTests<double, 2>();
    callSubTests<float, 3>();
    callSubTests<long double, 3>();
}
