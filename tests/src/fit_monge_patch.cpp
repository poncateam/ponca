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

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    // paraboloid parameters
    Scalar a = Eigen::internal::random<Scalar>(-2,2);
    Scalar b = Eigen::internal::random<Scalar>(-2,2);
    VectorType direction (0,0,1);

    std::cout << "plot(" << a << "x^2 + " << b << "y^2)" << std::endl;

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
//        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
//                                                     direction,
//                                                     width,
//                                                     _bAddPositionNoise,
//                                                     _bAddNormalNoise,
//                                                     _bUnoriented);

        vectorPoints[i] = getPointOnParaboloid<DataPoint>(a, b, analysisScale, _bAddPositionNoise);
    }

    epsilon = testEpsilon<Scalar>();
    if ( _bAddPositionNoise) // relax a bit the testing threshold
      epsilon = Scalar(0.02*MAX_NOISE);
    // Test for each point if the fitted plane correspond to the theoretical plane
#ifdef DEBUG
#pragma omp parallel for
#endif

    // evaluate only for the point located at 0,0, which is convenient as:
    //  - the normal vector is necessary 0,0,1
    //  - the tangent vectors are supposed to be 1,0,0 and 0,1,0
    //  - the fitted quadric should have the same parameters than the generate one
    const auto queryPos = VectorType::Zero();

    Fit fit;
    fit.setNeighborFilter({queryPos, analysisScale});
    fit.compute(vectorPoints);

    if( fit.isStable() ){
        std::cout << "plot("
                  << fit.h_uu() << "x^2 + "
                  << fit.h_vv() << "y^2 + "
                  << fit.h_uv() << "xy + "
                  << fit.h_u() << "x + "
                  << fit.h_v() << "y + "
                  << fit.h_c() <<")" << std::endl;

        // Check if the plane orientation is equal to the generation direction
        VERIFY(Scalar(1.) - std::abs(fit.plane().primitiveGradient(queryPos).dot(direction)) <= epsilon);

        // As the point cloud samples a plane, check if the quadric gradient is close to the plane's gradient
        VERIFY(Scalar(1.) -
               std::abs(fit.plane().primitiveGradient(queryPos).dot(fit.mongePatchPrimitive().primitiveGradient(queryPos))) <= epsilon);

        // Projecting to tangent plane and going back to world should not change the position
        VERIFY((fit.tangentPlaneToWorld(fit.worldToTangentPlane(queryPos)) - queryPos).norm() <= epsilon);

        if(!_bAddPositionNoise) {
          // Check if the query point is on the plane
          VERIFY(std::abs(fit.potential(queryPos)) <= epsilon);
          // check if we well have a plane
            auto first = fit.firstFundamentalForm();
            auto second = fit.secondFundamentalForm();
            auto w = fit.weingartenMap();
            auto kmean = fit.kMean();
            auto kmean2 = Scalar(0.5)*w.trace();
            auto gauss = fit.GaussianCurvature();
            auto gauss2 = w.determinant(); //return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_u()*u + h_v()*v + h_c()
          VERIFY(kmean <= epsilon);
          VERIFY(gauss <= epsilon);
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, MongePatchQuadraticFit> CovFitSmooth;
    typedef Basket<Point, WeightConstantFunc, MongePatchQuadraticFit> CovFitConstant;

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
