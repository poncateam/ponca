/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/gls_paraboloid_der.cpp
 \brief Test validity of GLS derivatives for a paraboloid
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc>
void testFunction()
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::QuaternionType QuaternionType;

    typedef typename Fit::VectorArray VectorArray;
    typedef typename Fit::ScalarArray ScalarArray;

    //generate sampled paraboloid
    int nbPoints = Eigen::internal::random<int>(1000, 10000);

    VectorType vCenter = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);
    VectorType vCoef = VectorType(Eigen::internal::random<Scalar>(-10,10), Eigen::internal::random<Scalar>(-10,10), 0);
    //vCoef.y() = vCoef.x();

    Scalar analysisScale = Scalar(.00000001 * std::min(fabs(vCoef.x()), fabs(vCoef.y())));


    Scalar rotationAngle = Eigen::internal::random<Scalar>(Scalar(0.), Scalar(2 * M_PI));
    VectorType vRotationAxis = VectorType::Random().normalized();
    QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
    qRotation = qRotation.normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar kappaEpsilon = 5e-1;

    vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnParaboloid<DataPoint>(vCenter, vCoef, qRotation, analysisScale, false);
    }

    Fit fit;
    fit.setWeightFunc(WeightFunc(analysisScale));
    VectorType vFittingPoint = VectorType(0, 0, 0);
    fit.init(vFittingPoint);

    for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
        it != vectorPoints.end();
        ++it)
    {
        fit.addNeighbor(*it);
    }

    fit.finalize();

    if(fit.isStable())
    {
        Scalar a = vCoef.x();
        Scalar b = vCoef.y();

        Scalar theoricTau = 0;
        VectorType theoricEta = VectorType(0, 0, 1);
        Scalar theoricKappa = (a + b) * Scalar(.5);

        Scalar computedTheoricKappa = getKappaMean<DataPoint>(vectorPoints, vFittingPoint, a, b, analysisScale);

        Scalar tau = fit.tau();
        VectorType eta = fit.eta();
        Scalar kappa = fit.kappa();

        VERIFY( Eigen::internal::isMuchSmallerThan(std::fabs(tau - theoricTau), 1., epsilon) );
        VERIFY( Eigen::internal::isMuchSmallerThan((theoricEta - eta).norm(), 1., epsilon ) );
        VERIFY( Eigen::internal::isMuchSmallerThan(std::fabs(computedTheoricKappa - kappa), 1., kappaEpsilon) );

        Scalar kappanorm = fit.kappa_normalized();
        Scalar taunorm = fit.tau_normalized();
        ScalarArray dkappa = fit.dkappa();

        Scalar kappa1 = fit.GLSk1();
        Scalar kappa2 = fit.GLSk2();
        Scalar meanKappaFromPricipalCurvatures = (kappa1 + kappa2) * Scalar(.5);

        //VERIFY( Eigen::internal::isApprox(meanKappaFromPricipalCurvatures, theoricKappa, kappaEpsilon) );

        Scalar gaussian = fit.GLSGaussianCurvature();
        Scalar theoricGaussian = a * b;

        //VERIFY( Eigen::internal::isApprox(gaussian, theoricGaussian, kappaEpsilon) );
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper> FitSmoothOriented;
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper> FitConstantOriented;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test paraboloid fitting..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
