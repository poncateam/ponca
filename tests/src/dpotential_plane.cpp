/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/dpotential_plane.cpp
    \brief Test validity of dPotential
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

    // Quick testing is requested for coverage
    int size = QUICK_TESTS ? 1 : int(vectorPoints.size());

    // Test for each point if the dPotential is correctly oriented
#pragma omp parallel for
    for(int i = 0; i < size; ++i)
    {
        Fit fit;
        fit.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            VectorType estimated = fit.dPotential().template tail<DataPoint::Dim>();

            Scalar norm = estimated.norm();
            Scalar absdot = std::abs(direction.dot(estimated/norm));

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

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, OrientedSphereSpaceDer> FitSmoothSpaceDer;
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, OrientedSphereSpaceDer> FitConstantSpaceDer;   
    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, OrientedSphereScaleSpaceDer> FitSmoothScaleSpaceDer;   
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, OrientedSphereScaleSpaceDer> FitConstantScaleSpaceDer;   
    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitSmoothSpaceMlsDer;
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitConstantSpaceMlsDer;   
    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitSmoothScaleSpaceMlsDer;   
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitConstantScaleSpaceMlsDer;   

    cout << "Testing with perfect plane (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceMlsDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceMlsDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceMlsDer>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceMlsDer>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceMlsDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceMlsDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceMlsDer>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceMlsDer>(true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test dPotential" << endl;

    //callSubTests<double, 2>();
    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
