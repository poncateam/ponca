/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/dnormal_plane.cpp
    \brief Test validity of dnormal
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

#include "Ponca/src/Fitting/basket.h"
#include "Ponca/src/Fitting/mean.h"
#include "Ponca/src/Fitting/weightFunc.h"
#include "Ponca/src/Fitting/weightKernel.h"

using namespace std;

template<typename DataPoint, typename FitA, typename FitB, typename WeightFunc> //, typename Fit, typename WeightFunction>
void compareFit(bool _bAddPositionNoise = false, bool /*_bAddNormalNoise */= false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);
    int nbPoints = Eigen::internal::random<int>(100, 1000);
    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Random() * centerScale;

    std::vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, false, false, false);
    }

    FitA fitA;
    FitB fitB;
    fitA.setWeightFunc(WeightFunc(VectorType::Zero(), analysisScale));
    // fitB.setWeightFunc(WeightKernel(1));
    // fitA.compute(vectorPoints);
    // fitB.compute(vectorPoints);
    // std::cout << fitA.barycenter() << std::endl; // Should be equal to center
    // std::cout << fitB.barycenter() << std::endl; // Should be equal to center
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef Ponca::DistWeightFunc<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncLocal;
    typedef Ponca::DistWeightFuncGlobal<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncGlobal;
    typedef Ponca::Basket<Point, WeightConstantFuncLocal, Ponca::MeanPosition> FitLocal;
    typedef Ponca::Basket<Point, WeightConstantFuncGlobal, Ponca::MeanPosition> FitGlobal;

    cout << "Testing the barycenter..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( compareFit<Point, FitLocal, FitGlobal, WeightConstantFuncLocal>( )));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test dnormal" << endl;

    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
