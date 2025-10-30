/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/normal_plane.cpp
    \brief Test validity of normal estimations
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
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
    if( _bAddPositionNoise ) // relax a bit the testing threshold
        epsilon = Scalar(0.002*MAX_NOISE);

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise);
    }

    // Test for each point if the normal is correct
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            VectorType estimated = fit.normal();

            Scalar norm = estimated.norm();
            Scalar absdot = std::abs(estimated.dot(direction));

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

    cout << "Testing with perfect plane (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceDer>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceDer>(true, true) ));

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

    //callSubTests<double, 2>();
    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
