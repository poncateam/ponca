/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
    \file test/Grenaille/scale_space_consistency.cpp
    \brief Test consistency of scale/space dPotential and dNormal
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename FitScaleDer, typename FitSpaceDer, typename FitScaleSpaceDer, typename NeighborFilter>
void checkConsistency(const vector<DataPoint> vectorPoints, typename DataPoint::Scalar analysisScale, typename DataPoint::Scalar epsilon)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;

    // input type check
    VERIFY(  FitScaleDer::isScaleDer()      && !FitScaleDer::isSpaceDer() );
    VERIFY( !FitSpaceDer::isScaleDer()      &&  FitSpaceDer::isSpaceDer() );
    VERIFY(  FitScaleSpaceDer::isScaleDer() &&  FitScaleSpaceDer::isSpaceDer() );

    const int nbPoints = vectorPoints.size();

    // test that dPotential and dNormal wrt scale, space and both are the same regardless of the derivation type
#pragma omp parallel for
    for(int i=0; i<nbPoints; ++i)
    {
        FitScaleDer fitScaleDer;
        fitScaleDer.setNeighborFilter(NeighborFilter(vectorPoints[i].pos(), analysisScale));
        fitScaleDer.compute(vectorPoints);

        FitSpaceDer fitSpaceDer;
        fitSpaceDer.setNeighborFilter(NeighborFilter(vectorPoints[i].pos(), analysisScale));
        fitSpaceDer.compute(vectorPoints);

        FitScaleSpaceDer fitScaleSpaceDer;
        fitScaleSpaceDer.setNeighborFilter(NeighborFilter(vectorPoints[i].pos(), analysisScale));
        fitScaleSpaceDer.compute(vectorPoints);

        VERIFY( fitScaleDer.isStable()==fitSpaceDer.isStable() && fitScaleDer.isStable()==fitScaleSpaceDer.isStable() );

        if(fitScaleDer.isStable() && fitSpaceDer.isStable() && fitScaleSpaceDer.isStable())
        {
            typename FitScaleDer::ScalarArray dPotentialScaleDer = fitScaleDer.dPotential();
            typename FitSpaceDer::ScalarArray dPotentialSpaceDer = fitSpaceDer.dPotential();
            typename FitScaleSpaceDer::ScalarArray dPotentialScaleSpaceDer = fitScaleSpaceDer.dPotential();

            typename FitSpaceDer::VectorArray dNormalSpaceDer = fitSpaceDer.dNormal();
            typename FitScaleSpaceDer::VectorArray dNormalScaleSpaceDer = fitScaleSpaceDer.dNormal();

            Scalar dt_potential1 = dPotentialScaleDer[0];
            Scalar dt_potential2 = dPotentialScaleSpaceDer[0];

            VectorType dx_potential1 = dPotentialSpaceDer.template tail<DataPoint::Dim>();
            VectorType dx_potential2 = dPotentialScaleSpaceDer.template tail<DataPoint::Dim>();

            MatrixType dx_normal1 = dNormalSpaceDer.template rightCols<DataPoint::Dim>();
            MatrixType dx_normal2 = dNormalScaleSpaceDer.template rightCols<DataPoint::Dim>();

            VERIFY( std::abs(dt_potential1-dt_potential2) < epsilon );
            VERIFY( (dx_potential1-dx_potential2).norm() < epsilon );

            for(int i=0; i<DataPoint::Dim; ++i)
            {
                VERIFY( (dx_normal1.col(i)-dx_normal2.col(i)).norm() < epsilon );
            }
        }
        else
        {
            cout << "Warning: unstable fit" << endl;
        }
    }
}



template<typename DataPoint, typename Fit1, typename Fit2, typename Fit3, typename NeighborFilter>
void testFunction(bool _bAddNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::QuaternionType QuaternionType;

    int nbPointsParaboloid = Eigen::internal::random<int>(100, 1000);
    int nbPointsPlane = Eigen::internal::random<int>(100, 1000);
    int nbPointsSphere = Eigen::internal::random<int>(100, 1000);
    int nbPointsRectangle = Eigen::internal::random<int>(100, 1000);

    VectorType coefParaboloid = 10 * VectorType(Eigen::internal::random<Scalar>(-1,1), Eigen::internal::random<Scalar>(-1,1), 0);
    Scalar analysisScaleParaboloid = Scalar(0.001);
    QuaternionType qParaboloid;

    Scalar widthPlane  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar heightPlane = widthPlane;
    Scalar analysisScalePlane = Scalar(50.) * std::sqrt( widthPlane * heightPlane / nbPointsPlane);
    Scalar centerScalePlane   = Eigen::internal::random<Scalar>(1,10000);
    VectorType centerPlane    = VectorType::Random() * centerScalePlane;
    VectorType directionPlane = VectorType::Random().normalized();

    Scalar radiusSphere = Eigen::internal::random<Scalar>(1,10);
    VectorType centerSphere = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);
    Scalar analysisScaleSphere = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radiusSphere * radiusSphere / nbPointsSphere);

    Scalar widthRectangle  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar heightRectangle = widthRectangle;
    Scalar analysisScaleRectangle = Scalar(20.) * std::sqrt( widthRectangle * heightRectangle / nbPointsRectangle);
    Scalar centerScaleRectangle   = Eigen::internal::random<Scalar>(1,10000);
    VectorType centerRectangle    = VectorType::Random() * centerScaleRectangle;
    VectorType directionRectangle = VectorType::Random().normalized();
    VectorType xAxisRectangle, yAxisRectangle;
    int i0 = -1, i1 = -1, i2 = -1;
    directionRectangle.minCoeff(&i0);
    i1 = (i0+1)%3;
    i2 = (i0+2)%3;
    xAxisRectangle[i0] = 0;
    xAxisRectangle[i1] = directionRectangle[i2];
    xAxisRectangle[i2] = -directionRectangle[i1];
    yAxisRectangle[i0] = directionRectangle[i1]*directionRectangle[i1] + directionRectangle[i2]*directionRectangle[i2];
    yAxisRectangle[i1] = -directionRectangle[i1]*directionRectangle[i0];
    yAxisRectangle[i2] = -directionRectangle[i2]*directionRectangle[i0];

    Scalar epsilon = testEpsilon<Scalar>();
    if( _bAddNoise ) // relax a bit the testing threshold
        epsilon = Scalar(0.002*MAX_NOISE);

    vector<DataPoint> vectorPointsParaboloid(nbPointsParaboloid);
    vector<DataPoint> vectorPointsPlane(nbPointsParaboloid);
    vector<DataPoint> vectorPointsSphere(nbPointsParaboloid);
    vector<DataPoint> vectorPointsRectangle(nbPointsRectangle);

    for(unsigned int i = 0; i < vectorPointsParaboloid.size(); ++i)
    {
        vectorPointsParaboloid[i] = getPointOnParaboloid<DataPoint>(VectorType::Zero(),
                                                                    coefParaboloid,
                                                                    qParaboloid,
                                                                    Scalar(1.2)*analysisScaleParaboloid,
                                                                    _bAddNoise);
    }
    for(unsigned int i = 0; i < vectorPointsPlane.size(); ++i)
    {
        vectorPointsPlane[i] = getPointOnPlane<DataPoint>(centerPlane,
                                                          directionPlane,
                                                          widthPlane,
                                                          _bAddNoise,
                                                          _bAddNoise);
    }
    for(unsigned int i = 0; i < vectorPointsParaboloid.size(); ++i)
    {
        vectorPointsSphere[i] = getPointOnSphere<DataPoint>(radiusSphere,
                                                            centerSphere,
                                                            _bAddNoise,
                                                            _bAddNoise);
    }
    for(unsigned int i = 0; i < vectorPointsRectangle.size(); ++i)
    {
        vectorPointsRectangle[i] = getPointOnRectangularPlane<DataPoint>(centerRectangle,
                                                                directionRectangle,
                                                                widthRectangle,
                                                                heightRectangle,
                                                                xAxisRectangle,
                                                                yAxisRectangle,
                                                                _bAddNoise);
        vectorPointsRectangle[i].normal() = directionRectangle;
    }

    checkConsistency<DataPoint, Fit1, Fit2, Fit3, NeighborFilter>(vectorPointsParaboloid, Scalar(1.2)*analysisScaleParaboloid, epsilon);
    checkConsistency<DataPoint, Fit1, Fit2, Fit3, NeighborFilter>(vectorPointsPlane, analysisScalePlane, epsilon);
    checkConsistency<DataPoint, Fit1, Fit2, Fit3, NeighborFilter>(vectorPointsSphere, analysisScaleSphere, epsilon);
    checkConsistency<DataPoint, Fit1, Fit2, Fit3, NeighborFilter>(vectorPointsRectangle, analysisScaleRectangle, epsilon);
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereScaleDer> FitConstantScaleDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereSpaceDer> FitConstantSpaceDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereScaleSpaceDer> FitConstantScaleSpaceDer;

    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereScaleDer, MlsSphereFitDer> FitConstantScaleMlsDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereSpaceDer, MlsSphereFitDer> FitConstantSpaceMlsDer;
    typedef Basket<Point,WeightConstantFunc,OrientedSphereFit,OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitConstantScaleSpaceMlsDer;

    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereScaleDer> FitSmoothScaleDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereSpaceDer> FitSmoothSpaceDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereScaleSpaceDer> FitSmoothScaleSpaceDer;

    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereScaleDer, MlsSphereFitDer> FitSmoothScaleMlsDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereSpaceDer, MlsSphereFitDer> FitSmoothSpaceMlsDer;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit,OrientedSphereScaleSpaceDer, MlsSphereFitDer> FitSmoothScaleSpaceMlsDer;

    cout << "Testing with perfect primitives..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleDer, FitConstantSpaceDer, FitConstantScaleSpaceDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantScaleMlsDer, FitConstantSpaceMlsDer, FitConstantScaleSpaceMlsDer, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleDer, FitSmoothSpaceDer, FitSmoothScaleSpaceDer, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothScaleMlsDer, FitSmoothSpaceMlsDer, FitSmoothScaleSpaceMlsDer, WeightSmoothFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Check consistency of scale/space methods" << endl;

    callSubTests<float, 3>();
    callSubTests<long double, 3>();
    callSubTests<double, 3>();
}
