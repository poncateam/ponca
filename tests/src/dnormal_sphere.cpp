/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/dnormal_sphere.cpp
    \brief Test validity of dnormal
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
    typedef typename DataPoint::MatrixType MatrixType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    Scalar curvature = Scalar(1.)/radius;

    VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = Eigen::internal::random<Scalar>(0.3, std::sqrt(2.f)) * radius;

    Scalar epsilon = testEpsilon<Scalar>();
    if( _bAddPositionNoise ) // relax a bit the testing threshold
        epsilon = Scalar(0.1);

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise);
    }

    // Test for each point if dNormal is correct
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setNeighborFilter(NeighborFilter(vectorPoints[i].pos(), analysisScale));
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            VectorType normal = (vectorPoints[i].pos()-center).normalized();

            MatrixType dN = fit.dNormal().template rightCols<DataPoint::Dim>();

            Eigen::SelfAdjointEigenSolver<MatrixType> solver(dN);

            VERIFY( std::abs(solver.eigenvalues()[0]) < epsilon );
            VERIFY( std::abs(solver.eigenvalues()[1]-curvature) < epsilon);
            VERIFY( std::abs(solver.eigenvalues()[2]-curvature) < epsilon);

            VERIFY( std::abs(1-std::abs(normal.dot(solver.eigenvectors().col(0)))) < epsilon );
            VERIFY( std::abs(normal.dot(solver.eigenvectors().col(1))) < epsilon );
            VERIFY( std::abs(normal.dot(solver.eigenvectors().col(2))) < epsilon );
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereSpaceDer> FitSmoothOrientedSpaceDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereSpaceDer> FitConstantOrientedSpaceDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer> FitSmoothOrientedScaleSpaceDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer> FitConstantOrientedScaleSpaceDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitSmoothOrientedSpaceMlsDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereSpaceDer, MlsSphereFitDer> FitConstantOrientedSpaceMlsDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitSmoothOrientedScaleSpaceMlsDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit, OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitConstantOrientedScaleSpaceMlsDer;

    cout << "Testing with perfect plane (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpaceMlsDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpaceMlsDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedScaleSpaceMlsDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedScaleSpaceMlsDer, WeightConstantFunc>() ));
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
