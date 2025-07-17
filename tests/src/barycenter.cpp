/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/barycenter.cpp
    \brief Test validity of Global and Local Weight Func
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

#include "Ponca/src/Fitting/basket.h"
#include "Ponca/src/Fitting/mean.h"
#include "Ponca/src/Fitting/weightFunc.h"
#include "Ponca/src/Fitting/weightKernel.h"

using namespace std;

template<typename DataPoint, typename FitA, typename FitB>
void compareFit(const bool _bAddPositionNoise = false, const bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    Scalar radius = Eigen::internal::random<int>(1, 10);
    int nbPoints = Eigen::internal::random<int>(100, 1000);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Zero();

    std::vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise, false);

    FitA fitA;
    FitB fitB;
    fitA.setNeighborFilter({center, radius});
    fitB.setNeighborFilter({center, radius});
    fitA.compute(vectorPoints);
    fitB.compute(vectorPoints);

    // std::cout << center.transpose() << std::endl;
    // std::cout << fitA.barycenter().transpose() << std::endl;
    // std::cout << fitB.barycenter().transpose() << std::endl;
    VERIFY(fitA.barycenter().isApprox(fitB.barycenter()));

    // TODO : Fix this test
    // Barycenter should also be somewhat close to the center
    // float eps =  0.1f; // Greater tolerance
    // VERIFY(center.isApprox(fitA.barycenter(), eps));
    // VERIFY(center.isApprox(fitB.barycenter(), eps));
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef Ponca::DistWeightFunc<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncLocal;
    typedef Ponca::DistWeightFuncGlobal<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncGlobal;
    typedef Ponca::Basket<Point, WeightConstantFuncLocal, Ponca::MeanPosition> FitConstantLocal;
    typedef Ponca::Basket<Point, WeightConstantFuncGlobal, Ponca::MeanPosition> FitConstantGlobal;

    typedef Ponca::DistWeightFunc<Point, Ponca::SmoothWeightKernel<Scalar> > WeightSmoothFuncLocal;
    typedef Ponca::DistWeightFuncGlobal<Point, Ponca::SmoothWeightKernel<Scalar> > WeightSmoothFuncGlobal;
    typedef Ponca::Basket<Point, WeightSmoothFuncLocal, Ponca::MeanPosition> FitSmoothLocal;
    typedef Ponca::Basket<Point, WeightSmoothFuncGlobal, Ponca::MeanPosition> FitSmoothGlobal;

    cout << "Testing the barycenter..." << endl;

    for(int i = 0; i < g_repeat; ++i) {
        CALL_SUBTEST(( compareFit<Point, FitConstantLocal, FitConstantGlobal>( )));
        CALL_SUBTEST(( compareFit<Point, FitSmoothLocal, FitSmoothGlobal>( )));
    }

    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test Global / Local Weight Func" << endl;

    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
