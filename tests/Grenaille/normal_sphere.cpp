/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/normal_sphere.cpp
    \brief Test validity of normal estimations
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius / nbPoints);

    Scalar epsilon = testEpsilon<Scalar>();
    //Scalar radisuEpsilon = epsilon * radius;

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise);
    }

    // Test for each point if the normal is correct
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints.cbegin(), vectorPoints.cend());

        if(fit.isStable())
        {
            VectorType estimated = fit.normal();
            VectorType theoritical = (vectorPoints[i].pos()-center).normalized();

            Scalar norm = estimated.norm();
            Scalar absdot = std::abs(estimated.dot(theoritical));

            VERIFY( std::abs(norm-Scalar(1.)) < epsilon );
            VERIFY( std::abs(absdot-Scalar(1.)) < epsilon );
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit> FitSmoothOriented;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit> FitConstantOriented;
    typedef Basket<Point,WeightSmoothFunc,UnorientedSphereFit> FitSmoothUnoriented;
    typedef Basket<Point,WeightConstantFunc,UnorientedSphereFit> FitConstantUnoriented;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitSmoothOrientedSpaceDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitConstantOrientedSpaceDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitSmoothOrientedScaleSpaceDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitConstantOrientedScaleSpaceDer;

    cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceDer, WeightSmoothFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceDer, WeightConstantFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceDer, WeightSmoothFunc>(true,true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceDer, WeightConstantFunc>(true,true) ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceDer, WeightConstantFunc>() ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test normal estimation" << endl;

    callSubTests<double, 2>();
    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
