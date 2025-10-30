/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/dpotential_sphere.cpp
    \brief Test validity of dPotential
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename NeighborFilter>
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

    // Test for each point if the dPotential is correctly oriented
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setNeighborFilter(NeighborFilter(vectorPoints[i].pos(), analysisScale));
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            VectorType estimated = fit.dPotential().template tail<DataPoint::Dim>();
            VectorType normal    = (center-vectorPoints[i].pos()).normalized();

            Scalar norm = estimated.norm();
            Scalar absdot = std::abs(normal.dot(estimated/norm));

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

    cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceMlsDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceMlsDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceMlsDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceMlsDer, WeightConstantFunc>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceDer, WeightSmoothFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceDer, WeightConstantFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceDer, WeightSmoothFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceDer, WeightConstantFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothSpaceMlsDer, WeightSmoothFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantSpaceMlsDer, WeightConstantFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleSpaceMlsDer, WeightSmoothFunc>(true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleSpaceMlsDer, WeightConstantFunc>(true, true) ));       
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
