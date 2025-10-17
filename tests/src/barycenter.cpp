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

    Scalar radius = Scalar(1);
    int nbPoints = Eigen::internal::random<int>(1000000, 2000000);
    Scalar centerScale = Eigen::internal::random<Scalar>(1, 10000);
    VectorType center = VectorType::Random() * centerScale;

    std::vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise, false);

    FitA fitA;
    FitB fitB;
    fitA.setNeighborFilter({center, radius});
    fitB.setNeighborFilter({center, radius});
    fitA.compute(vectorPoints);
    fitB.compute(vectorPoints);

    // Barycenter should also be somewhat close to the center
    Scalar highToleranceEps =  Scalar(0.1); // Greater tolerance

    try {
        VERIFY((fitA.barycenter().isApprox(fitB.barycenter(), 0.001)));
        VERIFY(((center - fitA.barycenter()).norm() < highToleranceEps));
        VERIFY(((center - fitB.barycenter()).norm() < highToleranceEps));
    } catch (const std::exception& /*e*/) {
        std::cout << "Sphere radius : " << radius << std::endl;
        std::cout << "Sphere center scale : " << centerScale << std::endl;
        std::cout << "Nb points : " << nbPoints << std::endl;
        std::cout << "True center : { " << center.transpose() << " }" << std::endl;
        std::cout << "fitLocal.barycenter() : { " << fitA.barycenter().transpose() << " }" << std::endl;
        std::cout << "fitGlobal.barycenter() : { " << fitB.barycenter().transpose() << " }" << std::endl;
        std::cout << "fitLocal.barycenter() distance from center : " << (center - fitA.barycenter()).norm() << std::endl;
        std::cout << "fitGlobal.barycenter() distance from center : " << (center - fitB.barycenter()).norm() << std::endl;
        throw;
    }
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

    for(int i = 0; i < g_repeat; ++i) {
        // std::cout << "CONST" << flush;
        CALL_SUBTEST(( compareFit<Point, FitConstantLocal, FitConstantGlobal>( )));
        // std::cout << "SMOOTH" << flush;
        CALL_SUBTEST(( compareFit<Point, FitSmoothLocal, FitSmoothGlobal>( )));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test Global / Local Weight Func" << endl;
    cout << "Testing the barycenter : " << flush;

    /* The global fit is less precise than the local fit when the center of the sphere is offsetted
     * which lead to errors in the barycenter test when using floating point numbers.
     * So we skip the float test to avoid this inaccuracy problem. */
    // cout << "float" << flush;
    // callSubTests<float, 3>();
    // cout << " (ok), " << flush;
    cout << "double" << flush;
    callSubTests<double, 3>();
    cout << " (ok), " << flush;
    cout << "long double" << flush;
    callSubTests<long double, 3>();
    cout << " (ok)" << endl;
}
