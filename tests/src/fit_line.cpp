
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 \file test/Grenaille/fit_line.cpp
 \brief Test validity of line fitting procedure(s)
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covarianceLineFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>
#include <iostream>

using namespace std;
using namespace Ponca;


template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bAddPositionNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled line
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    VectorType _direction = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar noiseMagnitude = 0.001;

    vector<DataPoint> vectorPoints(nbPoints);
    std::generate(vectorPoints.begin(), vectorPoints.end(), [&_direction,_bAddPositionNoise,noiseMagnitude]() {
        return DataPoint(_direction * Eigen::internal::random<Scalar>(0.1,2)
                + ( _bAddPositionNoise ? (VectorType::Random() * noiseMagnitude).eval() : VectorType::Zero() ));
    });

    epsilon = testEpsilon<Scalar>();

#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(1));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints.cbegin(), vectorPoints.cend());

        VERIFY( fit.isStable() );

        if(_bAddPositionNoise)
            VERIFY((fit.distance(vectorPoints[i].pos())) <= Scalar(4)*noiseMagnitude );
        else
            VERIFY((fit.distance(vectorPoints[i].pos())) <= epsilon);

        VERIFY( Scalar(1) - std::abs( fit.direction().dot(_direction) ) <= epsilon );
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPosition<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, CovarianceLineFit> LeastSquareFitSmooth;
    typedef Basket<Point, WeightConstantFunc, CovarianceLineFit> LeastSquareFitConstant;


    cout << "Testing with perfect line..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect line
        CALL_SUBTEST(( testFunction<Point, LeastSquareFitSmooth, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, LeastSquareFitConstant, WeightConstantFunc>() ));

    }
    cout << "Ok!" << endl;

     cout << "Testing with noise on position" << endl;
     for(int i = 0; i < g_repeat; ++i)
     {
         CALL_SUBTEST(( testFunction<Point, LeastSquareFitSmooth, WeightSmoothFunc>(true) ));
         CALL_SUBTEST(( testFunction<Point, LeastSquareFitConstant, WeightConstantFunc>(true) ));
     }
     cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test line fitting for different baskets..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
