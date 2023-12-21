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
#include <Ponca/src/Fitting/curvatureEstimation.h>
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
    //typedef typename Fit::WeightFunction WeightFunc;
    typedef typename Fit::DataPoint DataPoint;
    typedef typename TestFit::WeightFunction TestWeightFunc;
    typedef typename TestFit::DataPoint TestDataPoint;
    typedef typename TestDataPoint::Scalar TestScalar;
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    //typedef typename DataPoint::MatrixType MatrixType;

    //generate sampled paraboloid
    int nbPoints = Eigen::internal::random<int>(10000, 20000);

    // vCenter is ignored in getPointOnParaboloid
    VectorType vCenter = VectorType::Zero(); //::Random() * Eigen::internal::random<Scalar>(0, 1);
    VectorType vCoef = 100 * VectorType(Eigen::internal::random<Scalar>(-1,1), Eigen::internal::random<Scalar>(-1,1), 0);

    Scalar analysisScale {.001};// / std::max(std::abs(vCoef.x()), std::abs(vCoef.y()));
    vCenter *= analysisScale;

    Scalar rotationAngle = Eigen::internal::random<Scalar>(Scalar(0.), Scalar(2 * M_PI));
    VectorType vRotationAxis = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    Scalar approxEpsilon {0.1};

    vector<DataPoint> vectorPoints(nbPoints);
    vector<DataPoint> vectorPointsOrigin(nbPoints);
    vector<TestDataPoint> testVectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
      vectorPointsOrigin[i] = getPointOnParaboloid<DataPoint>(vCenter, vCoef, analysisScale*Scalar(1.2), false);
      // Add noise:
      // vectorPointsOrigin[i].pos() += VectorType::Random()*1e-6;
      //vectorPointsOrigin[i].normal() = (vectorPointsOrigin[i].normal() + VectorType::Random()*1e-6).normalized();
      vectorPoints[i].pos() = vectorPointsOrigin[i].pos() + vCenter;
      vectorPoints[i].normal() = vectorPointsOrigin[i].normal();

      testVectorPoints[i].pos()    = vectorPoints[i].pos().template cast<TestScalar>();
      testVectorPoints[i].normal() = vectorPoints[i].normal().template cast<TestScalar>();
    }

    VectorType theoricNormal = VectorType(0, 0, -1);

    TestFit fit;

    fit.setWeightFunc(TestWeightFunc(analysisScale));
    VectorType vFittingPoint = vCenter;
    fit.init(vFittingPoint.template cast<TestScalar>());
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
      typedef DistWeightFunc<RefPoint, SmoothWeightKernel<RefScalar> > RefWeightFunc;

      vector<RefPoint> refVectorPoints(nbPoints);
      for(unsigned int i = 0; i < vectorPoints.size(); ++i)
      {
        refVectorPoints[i].pos() = vectorPoints[i].pos().template cast<RefScalar>();
        refVectorPoints[i].normal() = vectorPoints[i].normal().template cast<RefScalar>();
      }

      // Centered fit:
      RefFit ref_fit;
      ref_fit.setWeightFunc(RefWeightFunc(analysisScale));
      VectorType vFittingPoint = vCenter;
      ref_fit.init(vFittingPoint.template cast<RefScalar>());
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
      Scalar h = Scalar(0.000001)*analysisScale;

      // Differentiation along scale, x, y, z:
      for(int k = 0; k<4; ++k)
      {
        RefFit f;
        f.setWeightFunc(RefWeightFunc(analysisScale));
        VectorType vFittingPoint = vCenter;
        if(k==0)
          f.setWeightFunc(RefWeightFunc(analysisScale+h));
        else
          vFittingPoint(k-1) += h;
        f.init(vFittingPoint.template cast<RefScalar>());
        f.compute(refVectorPoints);

        RefScalar flip_f   = (f.isSigned() || (f.primitiveGradient().dot(theoricNormal.template cast<RefScalar>()) > 0 )) ? RefScalar(1) : RefScalar(-1);

//         Scalar uc = f.m_uc;
//         typename RefFit::VectorType ul = f.m_ul;
//         typename RefFit::VectorType sumP = f.m_sumP;

        // Take into account centered basis:
//         if(k>0) uc        += -f.m_ul(k-1) * h + f.m_uq * h * h;
//         if(k>0) ul[k-1]   -= 2. * f.m_uq * h;
//         if(k>0) sumP[k-1] += f.getWeightSum() * h;

//         dSumP.col(k)  = ( f.m_cog      - ref_fit.m_cog ) / h;
//         dSumP.col(k)  = ( sumP      - ref_fit.m_sumP  ) / h;
//         dUc(k)        = ( uc        - ref_fit.m_uc    ) / h;
//         dUq(k)        = ( f.m_uq    - ref_fit.m_uq    ) / h;
//         dUl.col(k)    = ( ul        - ref_fit.m_ul    ) / h;
//         dTau(k)    = ( f.tau()   - ref_fit.tau()   ) / h;

        dPotential(k) = ( flip_f*f.potential() - flip_ref * ref_fit.potential()   ) / h;
        dN.col(k)     = ( flip_f*f.primitiveGradient().normalized() - flip_ref * ref_fit.primitiveGradient().normalized() ) / h;
//         dKappa(k)     = ( f.kappa() - ref_fit.kappa() ) / h;
      }

//       std::cout << "== Numerical differentiation ==\n";
//       std::cout << "dPotential: "  << dPotential << " == " << fit.dPotential() << " ; " << (dPotential.template cast<Scalar>()-flip_fit*fit.dPotential()).norm()/dPotential.norm() << "\n";
//       std::cout << "dN:\n"  << dN << "\n == \n" << flip_fit*fit.dNormal() << " ; " << (dN.template cast<Scalar>()-flip_fit*fit.dNormal()).norm()/dN.norm() << "\n";
//       std::cout << "eig(dN): " << Eigen::EigenSolver<typename DataPoint::MatrixType>(dN.template cast<Scalar>().template rightCols<3>()).eigenvalues().transpose() << "\n\n";

//       std::cout << "dKappa: " << dKappa << " == " << fit.dkappa() << " ; " << (dKappa.template cast<Scalar>()-fit.dkappa()).norm()/dKappa.norm() << "\n";
//       std::cout << "dUc: "  << dUc << " == " << fit.m_dUc << " ; " << (dUc.template cast<Scalar>()-fit.m_dUc).norm()/dUc.norm() << "\n";
//       std::cout << "dUq: "  << dUq << " == " << fit.m_dUq << " ; " << (dUq.template cast<Scalar>()-fit.m_dUq).norm()/dUq.norm() << "\n";
//       std::cout << "dTau: " << dTau << " == " << fit.dtau() << " ; " << (dTau.template cast<Scalar>()-fit.dtau()).norm()/dTau.norm() << "\n";
//       std::cout << "dUl:\n" << dUl << "\n == \n" << fit.m_dUl << " ; " << (dUl.template cast<Scalar>()-fit.m_dUl).norm()/dUl.norm() << "\n";
//       std::cout << "dSumP:\n" << dSumP << "\n == \n" << fit.m_dSumP << " ; " << (dSumP.template cast<Scalar>()-fit.m_dSumP).norm()/dSumP.norm() << "\n";
//       std::cout << "dSumP:\n" << dSumP << "\n == \n" << fit.m_dCog << " ; " << (dSumP.template cast<Scalar>()-fit.m_dCog).norm()/dSumP.norm() << "\n";

//       VERIFY( dUc.template cast<Scalar>().isApprox( fit.m_dUc, der_epsilon ) );
//       VERIFY( dUq.template cast<Scalar>().isApprox( fit.m_dUq, der_epsilon ) );
//       VERIFY( (dUl.template cast<Scalar>()-fit.m_dUl).norm() < fit.m_ul.norm() * der_epsilon );
      // VERIFY( dTau.template cast<Scalar>().isApprox( fit.dtau(), der_epsilon ) );
      VERIFY( dPotential.template cast<Scalar>().isApprox( flip_fit*fit.dPotential().template cast<Scalar>(), Scalar(der_epsilon) ) );

      VERIFY( dN.template cast<Scalar>().isApprox( flip_fit*fit.dNormal().template cast<Scalar>(), Scalar(der_epsilon) ) );
      //VERIFY( dKappa.template cast<Scalar>().isApprox( fit.dkappa(), der_epsilon ) );
    }

    VERIFY(fit.isStable());

    if(fit.isStable())
    {
      Scalar a = vCoef.x();
      Scalar b = vCoef.y();

      Scalar theoricPotential = 0;
//      Scalar theoricKmean        = (a + b) / Scalar(2.);
      Scalar theoricAverageKmean = getKappaMean<DataPoint>(vectorPointsOrigin, vFittingPoint, a, b, analysisScale);
      Scalar theoricK1 = (std::abs(a)<std::abs(b) ? b : a);
      Scalar theoricK2 = (std::abs(a)<std::abs(b) ? a : b);
      Scalar theoricGaussian = a * b;

      Scalar potential  = flip_fit * fit.potential();
      VectorType normal = flip_fit * fit.primitiveGradient().template cast<Scalar>();
//       Scalar kmean = fit.kappa();

      // principal curvatures k1,k2 are used here such that |k1| > |k2| (instead of kmin < kmax)
      Scalar kmin = fit.kmin();
      Scalar kmax = fit.kmax();
      Scalar kappa1 = flip_fit * (std::abs(kmin)<std::abs(kmax) ? kmax : kmin);
      Scalar kappa2 = flip_fit * (std::abs(kmin)<std::abs(kmax) ? kmin : kmax);
      Scalar kmeanFromK1K2 = (kappa1 + kappa2) * Scalar(.5);
      Scalar gaussian = fit.GaussianCurvature();

//       std::cout << "flip_fit : " << flip_fit << " , " << bool(isSigned || fit.primitiveGradient().dot(theoricNormal) > 0) << "\n";
//       std::cout << "potential : " << potential << "  \tref: " << theoricPotential << std::endl;
//       std::cout << "k1        : " << kappa1 << "  \tref: " << theoricK1 << std::endl;
//       std::cout << "k2        : " << kappa2 << "  \tref: " << theoricK2 << std::endl;
//       std::cout << "kmean     : " << /*kmean << ", " <<*/ kmeanFromK1K2 << "  \tref:" << theoricKmean << " , " << theoricAverageKmean << std::endl;
//       std::cout << "gaussian  : " << gaussian << "  \tref: " << theoricGaussian << std::endl;
//       std::cout << "normal    : " << normal.transpose() << "  \tref: " << theoricNormal.transpose() << std::endl;

      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(potential - theoricPotential), Scalar(1.), approxEpsilon) );
      VERIFY( Eigen::internal::isMuchSmallerThan((theoricNormal - normal).norm(), Scalar(1.), approxEpsilon ) );

      VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(theoricAverageKmean - kmeanFromK1K2), std::abs(theoricK1), approxEpsilon) );
//       VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(theoricAverageKmean - kmean), std::abs(theoricK1), approxEpsilon) );

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
            FitScaleSpaceDer, OrientedSphereDer, CurvatureEstimatorBase, NormalDerivativesCurvatureEstimator>;
    using RefFitSphereOriented = BasketDiff<
            Basket<RefPoint, RefWeightFunc, OrientedSphereFit>,
            FitScaleSpaceDer, OrientedSphereDer, CurvatureEstimatorBase, NormalDerivativesCurvatureEstimator>;
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
