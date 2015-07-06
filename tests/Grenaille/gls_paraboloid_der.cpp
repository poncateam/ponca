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

    //generate sampled paraboloid
    int nbPoints = Eigen::internal::random<int>(10000, 20000);

    VectorType vCenter = VectorType::Random() * Eigen::internal::random<Scalar>(1, 1000);
    VectorType vCoef = VectorType(Eigen::internal::random<Scalar>(-10,10), Eigen::internal::random<Scalar>(-10,10), 0);
    //vCoef.y() = vCoef.x();

    Scalar analysisScale = Scalar(.00001);// * std::min(std::abs(vCoef.x()), std::abs(vCoef.y()));

    Scalar rotationAngle = Eigen::internal::random<Scalar>(Scalar(0.), Scalar(2 * M_PI));
    VectorType vRotationAxis = VectorType::Random().normalized();
    QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
    qRotation = qRotation.normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar kappaEpsilon = 0.1;

    vector<DataPoint> vectorPoints(nbPoints);
    vector<DataPoint> vectorPointsOrigin(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
      vectorPointsOrigin[i] = getPointOnParaboloid<DataPoint>(vCenter, vCoef, qRotation, analysisScale, false);
      vectorPoints[i].pos() = qRotation * vectorPointsOrigin[i].pos();
      vectorPoints[i].normal() = qRotation * vectorPointsOrigin[i].normal();
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

    VERIFY(fit.isStable());
    
    if(fit.isStable())
    {
      Scalar a = vCoef.x();
      Scalar b = vCoef.y();

      Scalar theoricTau = 0;
      VectorType theoricNormal = qRotation * VectorType(0, 0, 1);
      Scalar theoricKmean        = (a + b) * Scalar(.5);
      Scalar theoricAverageKmean = getKappaMean<DataPoint>(vectorPointsOrigin, vFittingPoint, a, b, analysisScale);
      Scalar theoricK1 = std::abs(a)<std::abs(b) ? b : a;
      Scalar theoricK2 = std::abs(a)<std::abs(b) ? a : b;
      Scalar theoricGaussian = a * b;

      Scalar tau = fit.tau();
      VectorType normal = fit.eta();
      Scalar kmean = fit.kappa();

      Scalar kappa1 = fit.GLSk1();
      Scalar kappa2 = fit.GLSk2();
      Scalar kmeanFromK1K2 = (kappa1 + kappa2) * Scalar(.5);
      Scalar gaussian = fit.GLSGaussianCurvature();

//       std::cout << "k1        : " << kappa1 << "  \tref: " << theoricK1 << std::endl;
//       std::cout << "k2        : " << kappa2 << "  \tref: " << theoricK2 << std::endl;
//       std::cout << "kmean     : " << kmean << ", " << kmeanFromK1K2 << "  \tref:" << theoricKmean << " , " << theoricAverageKmean << std::endl;
//       std::cout << "gaussian  : " << gaussian << "  \tref: " << theoricGaussian << std::endl;
//       std::cout << "normal    : " << normal.transpose() << "  \tref: " << theoricNormal.transpose() << std::endl;
      
      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(tau - theoricTau), Scalar(1.), epsilon) );
      VERIFY( Eigen::internal::isMuchSmallerThan((theoricNormal - normal).norm(), Scalar(1.), epsilon ) );
      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(theoricAverageKmean - kmeanFromK1K2), std::abs(theoricK1), kappaEpsilon) );
      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(theoricAverageKmean - kmean), std::abs(theoricK1), kappaEpsilon) );
      
      if(std::abs(std::abs(theoricK1)-std::abs(theoricK2))>kappaEpsilon*std::abs(theoricK1))
      {
        VERIFY( Eigen::internal::isApprox(kappa1, theoricK1, kappaEpsilon) );
        VERIFY( std::abs(kappa2-theoricK2) < kappaEpsilon*std::abs(kappa1) );
      }
      else
      {
        VERIFY( Eigen::internal::isApprox(std::abs(kappa1), std::abs(theoricK1), kappaEpsilon) );
        VERIFY( std::abs(std::abs(kappa2)-std::abs(theoricK2)) < kappaEpsilon*std::abs(kappa1) );
      }

      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(gaussian-theoricGaussian), std::max(std::abs(theoricK1), std::abs(theoricGaussian)), 2*kappaEpsilon) );
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

    CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test paraboloid fitting..." << endl;

    for(int i = 0; i < g_repeat; ++i)
    {
      callSubTests<float, 3>();
      callSubTests<double, 3>();
      callSubTests<long double, 3>();
    }
}
