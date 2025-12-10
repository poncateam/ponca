/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/gls_paraboloid_der.cpp
 \brief Test validity of GLS derivatives for a paraboloid
 */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../common/testing.h"
#include "../common/scalar_precision_check.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/curvature.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;

template<typename Fit, typename RefFit, typename TestFit>
void testFunction()
{
    // Define related structure
    typedef typename Fit::DataPoint DataPoint;
    typedef typename TestFit::DataPoint TestDataPoint;
    typedef typename TestDataPoint::Scalar TestScalar;
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    //typedef typename DataPoint::MatrixType MatrixType;

    //generate sampled paraboloid
    int nbPoints = Eigen::internal::random<int>(10000, 20000);

    Scalar paraboloidA = Eigen::internal::random<Scalar>(-0.5,0.5);
    Scalar paraboloidB = Eigen::internal::random<Scalar>(-0.5,0.5);

    Scalar analysisScale = static_cast<Scalar>(.001);// / std::max(std::abs(vCoef.x()), std::abs(vCoef.y()));

    Scalar rotationAngle = Eigen::internal::random<Scalar>(Scalar(0.), Scalar(2 * M_PI));
    VectorType vRotationAxis = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar approxEpsilon = static_cast<Scalar>(.1);

    vector<DataPoint> vectorPoints(nbPoints);
    vector<DataPoint> vectorPointsOrigin(nbPoints);
    vector<TestDataPoint> testVectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
      vectorPointsOrigin[i] = getPointOnParaboloid<DataPoint>(paraboloidA, paraboloidB, 0, 0, 0, 0,analysisScale*Scalar(1.2), false);
      // Add noise:
      // vectorPointsOrigin[i].pos() += VectorType::Random()*1e-6;
      //vectorPointsOrigin[i].normal() = (vectorPointsOrigin[i].normal() + VectorType::Random()*1e-6).normalized();
      vectorPoints[i].pos() = vectorPointsOrigin[i].pos();
      vectorPoints[i].normal() = vectorPointsOrigin[i].normal();

      testVectorPoints[i].pos()    = vectorPoints[i].pos().template cast<TestScalar>();
      testVectorPoints[i].normal() = vectorPoints[i].normal().template cast<TestScalar>();
    }

    VectorType theoricNormal = VectorType(0, 0, -1);

    TestFit fit;
    const VectorType vFittingPoint = VectorType::Zero();
    fit.setNeighborFilter({vFittingPoint.template cast<TestScalar>(), analysisScale});
    fit.init();
    for(typename vector<TestDataPoint>::iterator it = testVectorPoints.begin();
        it != testVectorPoints.end();
        ++it)
    {
        fit.addNeighbor(*it);
    }
    fit.finalize();

    Scalar flip_fit = (fit.isSigned() || (fit.primitiveGradient().template cast<Scalar>().dot(theoricNormal) > 0 )) ? Scalar(1) : Scalar(-1);
    {
      // Check derivatives wrt numerical differentiation
      // Use long double for stable numerical differentiation
      typedef long double RefScalar;
      typedef PointPositionNormal<RefScalar, 3> RefPoint;

      vector<RefPoint> refVectorPoints(nbPoints);
      for(unsigned int i = 0; i < vectorPoints.size(); ++i)
      {
        refVectorPoints[i].pos() = vectorPoints[i].pos().template cast<RefScalar>();
        refVectorPoints[i].normal() = vectorPoints[i].normal().template cast<RefScalar>();
      }

      // Centered fit:
      RefFit ref_fit;
      ref_fit.setNeighborFilter({vFittingPoint.template cast<RefScalar>(), analysisScale});
      ref_fit.init();
      for(typename vector<RefPoint>::iterator it = refVectorPoints.begin();
          it != refVectorPoints.end();
          ++it)
      {
          ref_fit.addNeighbor(*it);
      }
      ref_fit.finalize();
      RefScalar flip_ref = (ref_fit.isSigned() || (ref_fit.primitiveGradient().dot(theoricNormal.template cast<RefScalar>()) > 0 )) ? RefScalar(1) : RefScalar(-1);

      RefScalar der_epsilon = epsilon*10;
      if(Eigen::internal::is_same<Scalar,float>::value)
        der_epsilon = epsilon*100;
      // FIXME check whether numerical accuracy can be improved with float
      typename RefFit::VectorArray /*dUl,*/ dN;
//      typename RefFit::VectorArray dSumP;
      typename RefFit::ScalarArray dPotential/*, dUc, dUq, dTau, dKappa*/;
      RefScalar h = RefScalar(0.000001)*RefScalar(analysisScale);

      // Differentiation along scale, x, y, z:
      for(int k = 0; k<4; ++k)
      {
        RefFit f;
        typename RefPoint::VectorType p = vFittingPoint.template cast<RefScalar>();
        auto scale = static_cast<RefScalar>(analysisScale);
        if(k==0)
          scale += h;
        else
          p(k-1) += h;
        f.setNeighborFilter({p, scale});
        f.compute(refVectorPoints);

        RefScalar flip_f   = (f.isSigned() || (f.primitiveGradient().dot(theoricNormal.template cast<RefScalar>()) > 0 )) ? RefScalar(1) : RefScalar(-1);
        dPotential(k) = ( flip_f*f.potential() - flip_ref * ref_fit.potential()   ) / h;
        dN.col(k)     = ( flip_f*f.primitiveGradient().normalized() - flip_ref * ref_fit.primitiveGradient().normalized() ) / h;
      }

//       std::cout << "== Numerical differentiation ==\n";
      VERIFY( dPotential.template cast<Scalar>().isApprox( flip_fit*fit.dPotential().template cast<Scalar>(), Scalar(der_epsilon) ) );
      VERIFY( dN.template cast<Scalar>().isApprox( flip_fit*fit.dNormal().template cast<Scalar>(), Scalar(der_epsilon) ) );
    }

    VERIFY(fit.isStable());

    if(fit.isStable())
    {

      Scalar theoricPotential = 0;

      Scalar theoricK1, theoricK2;
        {
            Scalar tK1 = Scalar(2) * paraboloidA;
            Scalar tK2 = Scalar(2) * paraboloidB;

            theoricK1 = flip_fit * (std::abs(tK1)<std::abs(tK2) ? tK2 : tK1);
            theoricK2 = flip_fit * (std::abs(tK1)<std::abs(tK2) ? tK1 : tK2);
        }
      Scalar theoricGaussian = theoricK1 * theoricK2;
      Scalar theoricKmean    = Scalar(0.5)*(theoricK1+theoricK2);

      Scalar potential  = flip_fit * fit.potential();
      VectorType normal = flip_fit * fit.primitiveGradient().template cast<Scalar>();

      // principal curvatures k1,k2 are used here such that |k1| > |k2| (instead of kmin < kmax)
      Scalar kmin = fit.kmin();
      Scalar kmax = fit.kmax();
      Scalar kappa1 = flip_fit * (std::abs(kmin)<std::abs(kmax) ? kmax : kmin);
      Scalar kappa2 = flip_fit * (std::abs(kmin)<std::abs(kmax) ? kmin : kmax);
      Scalar kmeanFromK1K2 = (kappa1 + kappa2) * Scalar(.5);
      Scalar gaussian = fit.GaussianCurvature();

      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(potential - theoricPotential), Scalar(1.), approxEpsilon) );
      VERIFY( Eigen::internal::isMuchSmallerThan((theoricNormal - normal).norm(), Scalar(1.), approxEpsilon ) );

      // The error in mean curvature estimation must be smaller in magnitude wrt the largest curvature
      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(theoricKmean - kmeanFromK1K2), std::abs(theoricK1), approxEpsilon));

      if(std::abs(std::abs(theoricK1)-std::abs(theoricK2))>approxEpsilon*std::abs(theoricK1))
      {
        // absolute curvatures are clearly different
        VERIFY( Eigen::internal::isApprox(kappa1, theoricK1, approxEpsilon) );
        VERIFY( std::abs(kappa2-theoricK2) < approxEpsilon*std::abs(kappa1) );
      }
      else
      {
        // absolute curvatures are close to each other and therefore their order should be ignored
        VERIFY( Eigen::internal::isApprox(std::abs(kappa1), std::abs(theoricK1), approxEpsilon) );
        VERIFY( std::abs(std::abs(kappa2)-std::abs(theoricK2)) < approxEpsilon*std::abs(kappa1) );
      }

      // The errors on k1 and k2 are expected to be of the same order, we thus compare the accuracy of k_gauss to k1^2
      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(gaussian-theoricGaussian), theoricK1*theoricK1, approxEpsilon) );
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef long double RefScalar;
    typedef PointPositionNormal<RefScalar, 3> RefPoint;
    typedef DistWeightFunc<RefPoint, SmoothWeightKernel<RefScalar> > RefWeightFunc;

    typedef ScalarPrecisionCheck<Scalar,RefScalar> TestScalar;
    TestScalar::check_enabled = false; // set it to true to track diverging computations
//    typedef PointPositionNormal<TestScalar, 3> TestPoint;
//    typedef DistWeightFunc<TestPoint, SmoothWeightKernel<TestScalar> > TestWeightFunc;

    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;

    using FitSphereOriented    = BasketDiff<
            Basket<Point, WeightSmoothFunc, OrientedSphereFit>,
            FitScaleSpaceDer, OrientedSphereDer, CurvatureEstimatorBaseDer, NormalDerivativeWeingartenEstimator, WeingartenCurvatureEstimatorDer>;
    using RefFitSphereOriented = BasketDiff<
            Basket<RefPoint, RefWeightFunc, OrientedSphereFit>,
            FitScaleSpaceDer, OrientedSphereDer, CurvatureEstimatorBaseDer, NormalDerivativeWeingartenEstimator, WeingartenCurvatureEstimatorDer>;
//    using TestFitSphereOriented = BasketDiff<Basket<TestPoint, TestWeightFunc, OrientedSphereFit>,
//            internal::FitScaleDer | internal::FitScaleDer, OrientedSphereDer, NormalDerivativesCurvatureEstimator>;

    CALL_SUBTEST(( testFunction<FitSphereOriented, RefFitSphereOriented, /*TestFitSphereOriented*/FitSphereOriented>() ));
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
