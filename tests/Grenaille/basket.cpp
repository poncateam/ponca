/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/Grenaille/basket.cpp
    \brief Test basket utility functions
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit>
void testFunction()
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename Fit::WFunctor WeightFunc;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 200);

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);

    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Random() * centerScale;

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, false, false, false);
    }

    // Test for each point if the fitted sphere correspond to the theorical sphere
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit1, fit2;

        // use compute function
        fit1.setWeightFunc(WeightFunc(analysisScale));
        fit1.init(vectorPoints[i].pos());
        fit1.compute(vectorPoints.cbegin(), vectorPoints.cend());

        // use addNeighbor
        fit2.setWeightFunc(WeightFunc(analysisScale));
        fit2.init(vectorPoints[i].pos());
        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit2.addNeighbor(*it);
        }

        fit2.finalize();

        // also test comparison operators
        VERIFY(fit1 == fit1);
        VERIFY(fit2 == fit2);
        VERIFY(fit1 == fit2);
        VERIFY(! (fit1 != fit1));
        VERIFY(! (fit1 != fit2));
        VERIFY(! (fit2 != fit2));
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    // We test only primitive functions and not the fitting procedure
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;

    typedef Basket<Point, WeightFunc, CovariancePlaneFit> Plane;
    typedef Basket<Point, WeightFunc, OrientedSphereFit> Sphere;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, Plane>() ));
        CALL_SUBTEST(( testFunction<Point, Sphere>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Basket functions in 3 dimensions..." << endl;
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
    cout << "Ok..." << endl;

    cout << "Test Basket functions in 4 dimensions..." << endl;
    callSubTests<float, 4>();
    callSubTests<double, 4>();
    callSubTests<long double, 4>();
    cout << "Ok..." << endl;

    cout << "Test Basket functions in 5 dimensions..." << endl;
    callSubTests<float, 5>();
    callSubTests<double, 5>();
    callSubTests<long double, 5>();
    cout << "Ok..." << endl;

    cout << "Test Basket functions in 42 dimensions..." << endl;
    callSubTests<float, 42>();
    callSubTests<double, 42>();
    callSubTests<long double, 42>();
    cout << "Ok..." << endl;
}
