/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
    \file test/Grenaille/gls_tau.cpp
    \brief Test validity GLS tau param
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(100, 1000);
    vector<DataPoint> vectorPoints(nbPoints);

    //Random plane parameters
    Scalar centerScale =  Eigen::internal::random<Scalar>(0, 10000);
    VectorType vCenter = VectorType::Random() * centerScale;
    VectorType vPlaneNormal = VectorType::Random().normalized();


    Scalar range = Eigen::internal::random<Scalar>(1, 10);
    Scalar epsilon = testEpsilon<Scalar>();
    if(_bAddPositionNoise || _bAddNormalNoise)
        epsilon *= 10;

    Scalar analysisScale = Scalar(100. * std::sqrt(4 * M_PI * range * range / nbPoints));

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        Scalar radius = Eigen::internal::random<Scalar>(-range, range);
        vectorPoints[i] = getPointOnPlane<DataPoint>(vCenter, vPlaneNormal, radius, _bAddPositionNoise, _bAddNormalNoise, _bUnoriented);
    }

    // Test for each point if the point moved from distance d correspond to tau
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        // Take a random distance to the plane, not too large to have few points in weightning analysis
        Scalar distanceToPlane = Eigen::internal::random<Scalar>(-range, range);
        VectorType vEvaluationPoint = vectorPoints[i].pos() + distanceToPlane * vPlaneNormal;

        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vEvaluationPoint);

        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }

        fit.finalize();

        if(fit.isStable())
        {
            Scalar fitTau = fit.tau();
            fitTau = fabs(fitTau);
            distanceToPlane = fabs(distanceToPlane);

            // Test Tau
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(distanceToPlane - fitTau), 1., epsilon) );
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam> FitSmoothOriented;
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam> FitConstantOriented;
    typedef Basket<Point, WeightSmoothFunc, UnorientedSphereFit, GLSParam> FitSmoothUnoriented;
    typedef Basket<Point, WeightConstantFunc, UnorientedSphereFit, GLSParam> FitConstantUnoriented;

    cout << "Testing with perfect plane (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>() ));
    }
    cout << "Ok..." << endl;

    /*cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
    CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(false, true, true) ));
    CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>(false, true, true) ));
    CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true, true, true) ));
    CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true, true, true) ));
    }
    cout << "Ok..." << endl;*/
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test GLS tau param coherance..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
