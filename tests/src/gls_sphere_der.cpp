/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/Grenaille/gls_sphere_der.cpp
    \brief Test validity of GLS derivatives for a sphere
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

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    typedef typename Fit::ScalarArray ScalarArray;

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

    // Test for each point if the Derivatives of kappa are equal to 0
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            ScalarArray dkappa = fit.dkappa();

            for(int i = 0; i < dkappa.size(); ++i)
            {
                VERIFY( Eigen::internal::isMuchSmallerThan(dkappa[i], Scalar(1.), epsilon) );
            }
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    using FitSmoothOriented = BasketDiff<
            Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam>,
            FitScaleSpaceDer, OrientedSphereDer, GLSDer>;


    cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
    }
    cout << "Ok!" << endl;

    /*cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(true, true) ));
    }
    cout << "Ok!" << endl;*/
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test sphere derivatives with GLSParam and OrientedSphereFit..." << endl;

    callSubTests<double, 2>();
    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
