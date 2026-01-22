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

    int nbPoints = Eigen::internal::random<int>(10000, 20000);
    // the radius is between 1 and 10
    const Scalar radius = Eigen::internal::random<Scalar>(1, 10);
    // the center is randomly located at maximum distance of 5 from the origin
    const VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 5);

    std::vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise, false);

    // the fittingScale must be large enough so that all the points included, otherwise it is not possible to compare
    // local and global fits, as local fits will discard all points whose distance is larger than radius.
    const Scalar fittingScale = center.norm() + Scalar(1.1)*radius; //radius is slightly increased to avoid approx errors

    FitA fitA;
    FitB fitB;
    fitA.setNeighborFilter({center, fittingScale}); // center will be ignored for global basis
    fitB.setNeighborFilter({center, fittingScale});
    fitA.compute(vectorPoints);
    fitB.compute(vectorPoints);

    // Barycenter should also be somewhat close to the center
    Scalar epsilon =  Scalar(0.01); // Greater tolerance

    try {
//        VERIFY((fitA.barycenter().isApprox(fitB.barycenter(), 0.001)));
        VERIFY(((fitA.barycenter() - fitB.barycenter()).norm() < epsilon));
        // here we have imprecision as the sphere is poorly sampled, so the barycenter can be quite far from sphere
        // center. To be conservative, we consider that it cannot be outside a sphere of radius 10% smaller than the
        // theoretical one, but colocated.
        VERIFY(((center - fitA.barycenter()).norm() < radius*Scalar(0.1)));
    } catch (const std::exception& /*e*/) {
        std::cout << "Sphere radius : " << radius << std::endl;
        std::cout << "Fit scale : " << fittingScale << std::endl;
        std::cout << "Nb points : " << nbPoints << std::endl;
        std::cout << "True center : { " << center.transpose() << " }" << std::endl;
        std::cout << "fitA.barycenter() : { " << fitA.barycenter().transpose() << " }" << std::endl;
        std::cout << "fitB.barycenter() : { " << fitB.barycenter().transpose() << " }" << std::endl;
        std::cout << "(fitA-fitB).barycenter() : { " << (fitA.barycenter()-fitB.barycenter()).transpose() << " }" << std::endl;
        std::cout << "(fitA-fitB).barycenter().norm() : { " << (fitA.barycenter()-fitB.barycenter()).norm() << " }" << std::endl;
        std::cout << "fitA.barycenter() distance from center : " << (center - fitA.barycenter()).norm() << std::endl;
        std::cout << "fitB.barycenter() distance from center : " << (center - fitB.barycenter()).norm() << std::endl;
        throw;
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef Ponca::PointPositionNormal<Scalar, Dim> Point;

    typedef Ponca::DistWeightFunc<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncLocal;
    typedef Ponca::NoWeightFuncGlobal<Point> NoWeightFuncGlobal;
    typedef Ponca::NoWeightFunc<Point> NoWeightFunc;

    typedef Ponca::Basket<Point, WeightConstantFuncLocal, Ponca::MeanPosition> FitConstantLocal;
    typedef Ponca::Basket<Point, NoWeightFunc, Ponca::MeanPosition> FitNoWeightLocal;
    typedef Ponca::Basket<Point, NoWeightFuncGlobal, Ponca::MeanPosition> FitNoWeightGlobal;

    typedef Ponca::DistWeightFunc<Point, Ponca::SmoothWeightKernel<Scalar> > WeightSmoothFuncLocal;
    typedef Ponca::Basket<Point, WeightSmoothFuncLocal, Ponca::MeanPosition> FitSmoothLocal;

    for(int i = 0; i < g_repeat; ++i) {
        // std::cout << "CONST" << flush;
        // check if, in local basis, no weight is equivalent to constant weight
        CALL_SUBTEST(( compareFit<Point, FitConstantLocal, FitNoWeightLocal>( )));
        // check that no weight is stable in both local and global frames
        CALL_SUBTEST(( compareFit<Point, FitNoWeightGlobal, FitNoWeightLocal>( )));
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
    cout << "float" << flush;
    callSubTests<float, 3>();
    cout << " (ok), " << flush;
    cout << "double" << flush;
    callSubTests<double, 3>();
    cout << " (ok), " << flush;
    cout << "long double" << flush;
    callSubTests<long double, 3>();
    cout << " (ok)" << endl;
}
