/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/gls_compare.cpp
 \brief Test coherance of gls compareTo
 */


#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/gls.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;

template<typename DataPoint, typename Fit, typename NeighborFilter>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius1 = Eigen::internal::random<Scalar>(1, 5);
    Scalar radius2 = Eigen::internal::random<Scalar>(10, 50);

    Scalar analysisScale = Scalar(Scalar(10.) * std::sqrt( Scalar(4.) * Scalar(M_PI) * radius2 * radius2 / nbPoints));

    VectorType center = VectorType::Zero();

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> sphere1(nbPoints);
    vector<DataPoint> sphere2(nbPoints);

    for(int i = 0; i < nbPoints; ++i)
    {
        sphere1[i] = getPointOnSphere<DataPoint>(radius1, center, _bAddPositionNoise, _bAddNormalNoise, _bUnoriented);
        sphere2[i] = getPointOnSphere<DataPoint>(radius2, center, _bAddPositionNoise, _bAddNormalNoise, _bUnoriented);
    }

    // Test for each point if the fitted sphere correspond to the theorical sphere
#pragma omp parallel for
    for(int i = 0; i < nbPoints - 1; ++i)
    {
        Fit fit1, fit2, fit3;
        fit1.setNeighborFilter(NeighborFilter(sphere1[i].pos(), analysisScale));
        fit2.setNeighborFilter(NeighborFilter(sphere1[i+1].pos(), analysisScale));
        fit3.setNeighborFilter(NeighborFilter(sphere2[i].pos(), analysisScale));

        fit1.compute(sphere1);
        fit2.compute(sphere1);
        fit3.compute(sphere2);

        if(fit1.isStable() && fit2.isStable() && fit3.isStable())
        {
            Scalar value1 = fit1.compareTo(fit2);
            Scalar value2 = fit1.compareTo(fit3);

            VERIFY( Eigen::internal::isMuchSmallerThan(value1, Scalar(1.), epsilon) );
            VERIFY( value2 > epsilon );
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam> FitSmoothOriented;

    cout << "Testing with perfect spheres (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect sphere
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(false, true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test GLS compareTo()..." << endl;

    callSubTests<float, 2>();
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 2>();
    callSubTests<long double, 3>();
}
