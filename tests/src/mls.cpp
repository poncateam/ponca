/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/mls.cpp
    \brief Test basket utility functions
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include "../split_test_helper.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>

#include <vector>

using namespace std;
using namespace Ponca;

template<typename DataPoint, typename Fit>
void testFunction() {
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    // Generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius / nbPoints);

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i) {
        // Add noise to the point cloud
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, true, true);
    }

    // Test for each point if the normal is correct
#pragma omp parallel for
    for (int i = 0; i < int(vectorPoints.size()); ++i)
    {
        VectorType pos = vectorPoints[i].pos();
        Fit fit;
        fit.setWeightFunc({pos, analysisScale});
        fit.compute(vectorPoints);

        Fit fitMLS;
        fitMLS.setWeightFunc({pos, analysisScale});
        fitMLS.computeMLS(vectorPoints, 1000);

        if(fit.isStable()) {
            // Tests the primitiveGradient
            const VectorType& estimated = fit.primitiveGradient(pos);
            const VectorType& estimatedMLS = fitMLS.primitiveGradient(pos);
            VectorType theoriticalNormal = (pos-center).normalized();

            Scalar absdot = std::abs(estimated.dot(theoriticalNormal));
            Scalar absdotMLS = std::abs(estimatedMLS.dot(theoriticalNormal));

            // Verify that mls gives better result than single projection when comparing with the theoretical values
            // By checking if absdotMLS is closer to 1 than asbdot (taking into account approximation error using the epsilon)
            // In other words : absdot >= absdotMLS >= 1 OR absdot <= absdotMLS <= 1
            VERIFY( (absdot + epsilon >= absdotMLS && absdotMLS >= 1 - epsilon) || (absdot - epsilon <= absdotMLS && absdotMLS <= 1 + epsilon));
        }
    }
}

template<typename Fit1, typename Fit2>
void isSamePlane(const Fit1& fit1, const Fit2& fit2) {
    const auto &plane1 = fit1.compactPlane();
    const auto &plane2 = fit2.compactPlane();

// Test we fit the same plane
    VERIFY(plane1.isApprox(plane2));
}

template<typename Fit1, typename Fit2>
void isSameSphere(const Fit1& fit1, const Fit2& fit2) {
    const auto &sphere1 = fit1.algebraicSphere();
    const auto &sphere2 = fit2.algebraicSphere();

    // Test we fit the same plane
    VERIFY(sphere1 == sphere2);
}

template<typename Fit1, typename Fit2>
void hasSamePlaneDerivatives(const Fit1& fit1, const Fit2& fit2) {
    // Get covariance
    const auto& dpot1 = fit1.covariancePlaneDer().dPotential();
    const auto& dpot2 = fit2.covariancePlaneDer().dPotential();
    const auto& dnor1 = fit1.covariancePlaneDer().dNormal();
    const auto& dnor2 = fit2.covariancePlaneDer().dNormal();

    // Test we compute the same derivatives
    VERIFY(dpot1.isApprox( dpot2 ));
    VERIFY(dnor1.isApprox( dnor2 ));
}

template<typename Scalar, int Dim>
void callSubTests()
{
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    //! [SpecializedPointType]

    // We test only primitive functions and not the fitting procedure
    //! [WeightFunction]
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar>> WeightFunc;
    //! [WeightFunction]
    typedef Basket<Point, WeightFunc, OrientedSphereFit>      Sphere;
    //! [PlaneFitType]
    typedef Basket<Point, WeightFunc, CovariancePlaneFit>      Plane;
    //! [PlaneFitType]

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef Basket<Point,WeightSmoothFunc,OrientedSphereFit> FitSmoothOriented;

    // //! [PlaneFitDerTypes]
    // using PlaneScaleDiff = BasketDiff<Plane, FitScaleDer, CovariancePlaneDer>;
    // using PlaneSpaceDiff = BasketDiff<Plane, FitSpaceDer, CovariancePlaneDer>;
    // using PlaneScaleSpaceDiff = BasketDiff<Plane, FitScaleSpaceDer, CovariancePlaneDer>;
    // //! [PlaneFitDerTypes]

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, Sphere>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Basket functions in 3 dimensions: float" << flush;
    CALL_SUBTEST_1((callSubTests<float, 3>()));
    cout << " (ok), double" << flush;
    CALL_SUBTEST_2((callSubTests<double, 3>()));
    cout << " (ok)" << flush;
    cout << ", long double" << flush;
    CALL_SUBTEST_3((callSubTests<long double, 3>()));
    cout << " (ok)" << flush;

    // cout << "Test Basket functions in 4 dimensions..." << endl;
    // callSubTests<float, 4>();
    // callSubTests<double, 4>();
    // callSubTests<long double, 4>();
}
