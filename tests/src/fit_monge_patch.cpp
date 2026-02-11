/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/fit_monge_patch.cpp
 \brief Test validity of monge patch fitting procedure(s)
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
//#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/mongePatch.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>
#include <fstream>

using namespace std;
using namespace Ponca;


template<typename DataPoint, typename Fit>
void testFunction(bool _bAddPositionNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);

    // we want a large scale to avoid fitting errors
    Scalar analysisScale = Scalar(4.) * width;//std::sqrt( width * height / nbPoints);

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    // paraboloid parameters
    // we need to use small magnitude for quadratic terms, otherwise the plane fitting is going to be wrong
    const Scalar standardMagnitude = 1;
    const Scalar smallMagnitude = 0.1;
    auto quadParams = Eigen::Matrix<Scalar, 6, 1>(
            Eigen::internal::random<Scalar>(-smallMagnitude,smallMagnitude),
            Eigen::internal::random<Scalar>(-smallMagnitude,smallMagnitude),
            Eigen::internal::random<Scalar>(-smallMagnitude,smallMagnitude),
            Eigen::internal::random<Scalar>(-standardMagnitude,standardMagnitude),
            Eigen::internal::random<Scalar>(-standardMagnitude,standardMagnitude),
            Eigen::internal::random<Scalar>(-standardMagnitude,standardMagnitude));

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnParaboloid<DataPoint>(quadParams, width, _bAddPositionNoise);
    }

    epsilon = 0.1*width;
    if ( _bAddPositionNoise) // relax a bit the testing threshold
      epsilon = std::max(Scalar(0.5*MAX_NOISE), epsilon*2);
    // Test for each point if the fitted plane correspond to the theoretical plane


    // evaluate only for the point located at 0,0
    const auto queryPos = VectorType::Zero();

    Fit fit;
    fit.setNeighborFilter({queryPos, analysisScale});
    fit.compute(vectorPoints);

    VERIFY( fit.isStable() );
    {
        // compute RMSE
        Scalar error {0};
        int nbTestSamples = vectorPoints.size()/5; // test on fewer samples to speed things up

        std::vector<VectorType> testPoints(nbTestSamples);

        for(unsigned int i = 0; i < nbTestSamples; ++i)
        {
            // get samples on points closer to the center, to avoid over-estimated error near the border
            const auto point = getPointOnParaboloid<DataPoint>(quadParams, width/2, false).pos();
            error += (point - fit.project(point)).norm();
            testPoints[i]=point;
        }
        error /= nbTestSamples;

        if(error > epsilon) {
            std::cerr << "Precision test failed (" << error << "). Dumping files\n";
#ifndef PONCA_COVERAGE_ENABLED
            {
                std::ofstream inFile, projFile;
                inFile.open ("input.xyz");
                projFile.open ("projected.xyz");
                for(unsigned int i = 0; i < testPoints.size(); ++i)
                {
                    inFile << testPoints[i].transpose() << endl;
                    projFile << fit.project(testPoints[i]).transpose() << std::endl;
                }
                inFile.close();
                projFile.close();
            }
#endif
            VERIFY(error <= epsilon);
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPosition<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, MongePatchQuadraticFit> CovFitSmooth;
    typedef Basket<Point, WeightConstantFunc, MongePatchQuadraticFit> CovFitConstant;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with plane noise on position" << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth>(true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant>(true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Monge Patch fitting for different baskets..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
