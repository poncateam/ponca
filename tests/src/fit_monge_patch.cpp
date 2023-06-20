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
#include <Ponca/src/Fitting/mean.h>
#include <Ponca/src/Fitting/plane.h>
#include <Ponca/src/Fitting/localFrame.h>
#include <Ponca/src/Fitting/covarianceFit.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/mongePatch.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;


template<typename DataPoint, typename Fit>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar height = width;

    Scalar analysisScale = Scalar(15.) * std::sqrt( width * height / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1,10000);
    VectorType center    = VectorType::Random() * centerScale;

    VectorType direction = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise,
                                                     _bUnoriented);
    }

    epsilon = testEpsilon<Scalar>();
    if ( _bAddPositionNoise) // relax a bit the testing threshold
      epsilon = Scalar(0.02*MAX_NOISE);
    // Test for each point if the fitted plane correspond to the theoretical plane

    // Quick testing is requested for coverage
    int size = QUICK_TESTS ? 1 : int(vectorPoints.size());

#ifdef DEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < size; ++i)
    {
        const auto& queryPos = vectorPoints[i].pos();

        Fit fit;
        fit.setNeighborFilter({queryPos, analysisScale});
        fit.compute(vectorPoints);

        if( fit.isStable() ){

            // Check if the plane orientation is equal to the generation direction
            VERIFY(Scalar(1.) - std::abs(fit.primitiveGradient(queryPos).dot(direction)) <= epsilon);

            // Projecting to tangent plane and going back to world should not change the position
            VERIFY((fit.localFrameToWorld(fit.worldToLocalFrame(queryPos)) - queryPos).norm() <= epsilon);

            if(!_bAddPositionNoise) {
              // Check if the query point is on the plane
              VERIFY(std::abs(fit.potential(queryPos)) <= epsilon);
              // check if we well have a plane
              VERIFY(fit.kMean() <= epsilon);
              VERIFY(fit.GaussianCurvature() <= epsilon);
            }
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, CovariancePlaneFit, MongePatch> CovFitSmooth;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit, MongePatch> CovFitConstant;

    // \todo Add these tests when MeanPlaneFit PROVIDES_TANGENT_PLANE_BASIS
//    typedef Basket<Point, WeightSmoothFunc, MeanPlaneFit, MongePatch> MeanFitSmooth;
//    typedef Basket<Point, WeightConstantFunc, MeanPlaneFit, MongePatch> MeanFitConstant;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect plane
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant>() ));
//        CALL_SUBTEST(( testFunction<Point, MeanFitSmooth>() ));
//        CALL_SUBTEST(( testFunction<Point, MeanFitConstant>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with plane noise on position" << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant>(false, true, true) ));
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
