
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

#define MULTIPASS_LINE_FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/leastSquareLine.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>
#include <iostream>

using namespace std;
using namespace Ponca;


template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled line
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    VectorType _direction = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = {i * _direction, {0,0,0}}; 
    }

    epsilon = testEpsilon<Scalar>();
    // if ( _bAddPositionNoise) // relax a bit the testing threshold
    //   epsilon = Scalar(0.01*MAX_NOISE);
    // Test for each point if the fitted line correspond to the theoretical line
#ifdef DEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {

        Fit fit;
        fit.setWeightFunc(WeightFunc(0.1));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints.cbegin(), vectorPoints.cend());

        if( fit.isStable() ){

            // Check if the query point is on the line
               VERIFY(fit.distance(vectorPoints[i].pos()) <= epsilon);

            // Check if the line orientation is equal to the generation direction
            //VERIFY(Scalar(1.) - std::abs(fit.direction(vectorPoints[i].pos()).dot(_direction)) <= epsilon);
            // if(!_bAddPositionNoise)

        }
        else {
            VERIFY(MULTIPASS_LINE_FITTING_FAILED);
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, LeastSquareLine> LeastSquareFitSmooth;
    typedef Basket<Point, WeightConstantFunc, LeastSquareLine> LeastSquareFitConstant;


    cout << "Testing with perfect line..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect line
        CALL_SUBTEST(( testFunction<Point, LeastSquareFitSmooth, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, LeastSquareFitConstant, WeightConstantFunc>() ));

    }
    cout << "Ok!" << endl;

    // cout << "Testing with noise on position" << endl;
    // for(int i = 0; i < g_repeat; ++i)
    // {
    //     CALL_SUBTEST(( testFunction<Point, LeastSquareFitSmooth, WeightSmoothFunc>(false, true, true) ));
    //     CALL_SUBTEST(( testFunction<Point, LeastSquareFitConstant, WeightConstantFunc>(false, true, true) ));
    // }
    // cout << "Ok!" << endl;
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
