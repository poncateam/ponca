/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/mls.cpp
    \brief Test basket utility functions
 */

#include "Ponca/src/Fitting/mls.h"

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

template<typename DataPoint>
typename DataPoint::Scalar generateData(KdTree<DataPoint>& tree)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
#ifdef NDEBUG
    int nbPoints = Eigen::internal::random<int>(500, 1000);
#else
    int nbPoints = Eigen::internal::random<int>(100, 200);
#endif

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);

    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Random() * centerScale;

    vector<DataPoint> vectorPoints(nbPoints);

#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, false, false, false);
    }

    tree.clear();
    tree.build(vectorPoints);

    return analysisScale;
}
template<typename Fit1, typename Fit2, typename Functor>
void testMLSIsSame(const KdTree<typename Fit1::DataPoint>& tree,
                typename Fit1::Scalar analysisScale,
                Functor f)
{
    static_assert(std::is_same<typename Fit1::DataPoint, typename Fit2::DataPoint>::value, "Both Fit should use the same point type");
    static_assert(std::is_same<typename Fit1::WFunctor, typename Fit2::WFunctor>::value, "Both Fit should use the same WFunctor");

    // Define related structure
    typedef typename Fit1::Scalar     Scalar;
    typedef typename Fit1::VectorType VectorType;
    typedef typename Fit1::WFunctor   WeightFunc;
    const auto& vectorPoints = tree.points();

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        Fit1 fit1;
        Fit2 fit2;

        auto neighborhoodRange = tree.range_neighbors(vectorPoints[i].pos(), analysisScale);

        // use compute function
        fit1.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        fit1.computeWithIds( neighborhoodRange, vectorPoints );

        fit2.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        fit2.setIterMLS(3);
        fit2.computeWithIdsMLS( neighborhoodRange, vectorPoints );

        f(fit1, fit2);
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
    typedef Basket<Point, WeightFunc, OrientedSphereFit, MLS> SphereMLS;
    //! [PlaneFitType]
    typedef Basket<Point, WeightFunc, CovariancePlaneFit>      Plane;
    typedef Basket<Point, WeightFunc, CovariancePlaneFit, MLS> PlaneMLS;
    //! [PlaneFitType]

    //
    // //! [PlaneFitDerTypes]
    // using PlaneScaleDiff = BasketDiff<Plane, FitScaleDer, CovariancePlaneDer>;
    // using PlaneSpaceDiff = BasketDiff<Plane, FitSpaceDer, CovariancePlaneDer>;
    // using PlaneScaleSpaceDiff = BasketDiff<Plane, FitScaleSpaceDer, CovariancePlaneDer>;
    // //! [PlaneFitDerTypes]

    KdTreeDense<Point> tree;
    Scalar scale = generateData(tree);

    for(int i = 0; i < g_repeat; ++i)
    {

        //  Plane diffs
        // CALL_SUBTEST((testBasicFunctionalities<PlaneScaleDiff>(tree, scale) ));
        // CALL_SUBTEST((testBasicFunctionalities<PlaneSpaceDiff>(tree, scale) ));
        // CALL_SUBTEST((testBasicFunctionalities<PlaneScaleSpaceDiff>(tree, scale) ));

        // // Check that we get the same Sphere, whatever the extensions
        auto checkIsSameSphere = [](const auto&f1, const auto&f2){isSameSphere(f1,f2);};
        CALL_SUBTEST((testMLSIsSame<Sphere, SphereMLS>(tree, scale, checkIsSameSphere) ));
        //
        // // Check that we get the same Plane, whatever the extensions
        auto checkIsSamePlane = [](const auto&f1, const auto&f2){isSamePlane(f1,f2);};
        CALL_SUBTEST((testMLSIsSame<Plane, PlaneMLS>(tree, scale, checkIsSamePlane) ));
        //
        // auto checkIsSamePlaneDerivative = [](const auto&f1, const auto&f2){
        //     isSamePlane(f1,f2);
        //     hasSamePlaneDerivatives(f1, f2);
        // };
        // // // Check that we get the same Plane derivative, whatever if we have one or more primitive fitted.
        // // int fff = 0;
        // CALL_SUBTEST((testMLSIsSame<Plane, PlaneMLS>(tree, scale, checkIsSamePlaneDerivative) ));
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