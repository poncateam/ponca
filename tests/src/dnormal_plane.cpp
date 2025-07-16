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

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool _bAddPositionNoise = false, bool /*_bAddNormalNoise */= false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar height = width;

    Scalar analysisScale = Scalar(20.) * std::sqrt( width * height / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1,10000);
    VectorType center    = VectorType::Random() * centerScale;
    VectorType direction = VectorType::Random().normalized();
    VectorType xAxis, yAxis;
    int i0 = -1, i1 = -1, i2 = -1;
    direction.minCoeff(&i0);
    i1 = (i0+1)%3;
    i2 = (i0+2)%3;
    xAxis[i0] = 0;
    xAxis[i1] = direction[i2];
    xAxis[i2] = -direction[i1];
    yAxis[i0] = direction[i1]*direction[i1] + direction[i2]*direction[i2];
    yAxis[i1] = -direction[i1]*direction[i0];
    yAxis[i2] = -direction[i2]*direction[i0];

    Scalar epsilon = testEpsilon<Scalar>();
    if( _bAddPositionNoise ) // relax a bit the testing threshold
        epsilon = Scalar(.1);

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnRectangularPlane<DataPoint>(center,
                                                                direction,
                                                                width,
                                                                height,
                                                                xAxis,
                                                                yAxis,
                                                                _bAddPositionNoise);
        vectorPoints[i].normal() = direction;
    }

    // Test for each point if dNormal is null
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            MatrixType dN = fit.dNormal().template rightCols<DataPoint::Dim>();
            Scalar diff = dN.array().abs().maxCoeff();

            VERIFY( diff < epsilon );
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
