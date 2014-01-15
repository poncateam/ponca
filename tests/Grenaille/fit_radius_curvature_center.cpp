/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/fit_radius_curvature_center.cpp
 \brief Test validity of algebraic sphere procedure and GLS kappa
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

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);

    Scalar analysisScale = 10.f * std::sqrt( 4.f * M_PI * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Random() * centerScale;

    Scalar epsilon = testEpsilon<Scalar>();
    // epsilon is relative to the radius size
    Scalar radiusEpsilon = epsilon * radius;


    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise, _bUnoriented);
    }

    // Test for each point if the fitted sphere correspond to the theorical sphere
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());

        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }


        fit.finalize();

        if(fit.isStable())
        {
            Scalar fitRadiusKappa = Scalar(fabs(Scalar(1.) / fit.kappa()));
            Scalar fitRadiusAlgebraic = fit.radius();
            VectorType fitCenter = fit.center();

            Scalar radiusMax = radius * MAX_NOISE;
            Scalar radiusMin = radius * MIN_NOISE;

            // Test procedure
            VERIFY( (fitCenter - center).norm() < (radiusMax - radius) + radiusEpsilon );
            VERIFY( (fitRadiusAlgebraic > radiusMin - radiusEpsilon) && (fitRadiusAlgebraic < radiusMax + radiusEpsilon) );
            // Test reparametrization
            VERIFY( (fitRadiusKappa > radiusMin - radiusEpsilon) && (fitRadiusKappa < radiusMax + radiusEpsilon) );
            //Test coherance
            VERIFY( Eigen::internal::isMuchSmallerThan(std::fabs(fitRadiusAlgebraic - fitRadiusKappa), 1., epsilon) );

            //Test on eta
            if(!_bAddPositionNoise && !_bAddNormalNoise)
            {
                //sometimes eta can be reversed
                VectorType fitEta = fit.eta().normalized().array().abs();
                VectorType theoricEta = vectorPoints[i].normal().array().abs();

                VERIFY( Eigen::internal::isMuchSmallerThan((fitEta - theoricEta).norm(), 1., epsilon)  );
            }
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPosistionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam> FitSmoothOriented;
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam> FitConstantOriented;
    typedef Basket<Point, WeightSmoothFunc, UnorientedSphereFit, GLSParam> FitSmoothUnoriented;
    typedef Basket<Point, WeightConstantFunc, UnorientedSphereFit, GLSParam> FitConstantUnoriented;

    cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect sphere
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true) ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true, true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test sphere fitting (radius / center) and GLS curvature for different baskets..." << endl;

    callSubTests<float, 2>();
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 2>();
    callSubTests<long double, 3>();
}
