/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/basket.cpp
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

template<typename Fit>
void testBasicFunctionalities(const KdTree<typename Fit::DataPoint>& tree, typename Fit::Scalar analysisScale)
{
    using DataPoint = typename Fit::DataPoint;

    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename Fit::NeighborFilter NeighborFilter;

    const auto& vectorPoints = tree.points();

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        const auto &fitInitPos = vectorPoints[i].pos();

        // use addNeighbor
        //! [Fit Manual Traversal]
        Fit fit1;
        fit1.setNeighborFilter({fitInitPos, analysisScale});
        fit1.init();
        for(auto it = vectorPoints.begin(); it != vectorPoints.end(); ++it)
           fit1.addNeighbor(*it);
        fit1.finalize();
        //! [Fit Manual Traversal]

        // use compute function
        //! [Fit Compute]
        Fit fit2;
        fit2.setNeighborFilter({fitInitPos, analysisScale});
        fit2.compute(vectorPoints);
        //! [Fit Compute]

        // also test comparison operators
        VERIFY(fit1 == fit1);
        VERIFY(fit2 == fit2);
        VERIFY(fit1 == fit2);
        VERIFY(fit2 == fit1);
        VERIFY(! (fit1 != fit1));
        VERIFY(! (fit1 != fit2));
        VERIFY(! (fit2 != fit2));

        Fit fit3;
        fit3.setNeighborFilter({fitInitPos, analysisScale});
        // Sort fit1
        std::list<int> neighbors3;
        for (int iNeighbor : tree.range_neighbors(fitInitPos, analysisScale))
            neighbors3.push_back(iNeighbor);
        neighbors3.sort();
        // Compute the neighbors
        fit3.computeWithIds( neighbors3, vectorPoints );

        VERIFY(fit3 == fit3);
        VERIFY(fit1 == fit3);
        VERIFY(! (fit1 != fit3));
        VERIFY(fit3 == fit1);
    }
}

template<typename Fit1, typename Fit2, typename Functor>
void testIsSame(const KdTree<typename Fit1::DataPoint>& tree,
                typename Fit1::Scalar analysisScale,
                Functor f)
{
    static_assert(std::is_same<typename Fit1::DataPoint, typename Fit2::DataPoint>::value, "Both Fit should use the same point type");
    static_assert(std::is_same<typename Fit1::NeighborFilter, typename Fit2::NeighborFilter>::value, "Both Fit should use the same NeighborFilter");

    // Define related structure
    typedef typename Fit1::Scalar         Scalar;
    typedef typename Fit1::VectorType     VectorType;
    typedef typename Fit1::NeighborFilter NeighborFilter;
    const auto& vectorPoints = tree.points();

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        using Fit = Fit1; // create an alias to ease doc snippet generation
        //! [Fit computeWithIds]
        Fit fit3;
        fit3.setWeightFunc({vectorPoints[i].pos(), analysisScale});
        auto neighborhoodRange = tree.range_neighbors(vectorPoints[i].pos(), analysisScale);
        fit3.computeWithIds( neighborhoodRange, vectorPoints );
        //! [Fit computeWithIds]

        Fit2 fit2;
        fit2.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit2.computeWithIds( neighborhoodRange, vectorPoints );

        f(fit3, fit2);
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    //! [SpecializedPointType]

    // We test only primitive functions and not the fitting procedure
    //! [NeighborFilter]
    using NeighborFilter = DistWeightFunc<Point, SmoothWeightKernel<Scalar> >;
    //! [NeighborFilter]
    using Sphere     = Basket<Point, NeighborFilter, OrientedSphereFit>;
    //! [PlaneFitType]
    using TestPlane = Basket<Point, NeighborFilter, CovariancePlaneFit>;
    //! [PlaneFitType]
    //! [FitType]
    using Sphere = Basket<Point, NeighborFilter, OrientedSphereFit>;
    //! [FitType]
    //! [HybridType]
    // Create an hybrid structure fitting a plane and a sphere at the same time
    using Hybrid = Basket<Point, NeighborFilter,
                          AlgebraicSphere, Plane,                     // primitives
                          MeanNormal, MeanPosition,                   // shared computation
                          OrientedSphereFitImpl,                      // sphere fitting
                          CovarianceFitBase, CovariancePlaneFitImpl>; // plane fitting
    //! [HybridType]

    //! [PlaneFitDerTypes]
    using PlaneScaleDiff = BasketDiff<TestPlane, FitScaleDer, CovariancePlaneDer>;
    using PlaneSpaceDiff = BasketDiff<TestPlane, FitSpaceDer, CovariancePlaneDer>;
    using PlaneScaleSpaceDiff = BasketDiff<TestPlane, FitScaleSpaceDer, CovariancePlaneDer>;
    //! [PlaneFitDerTypes]

    using HybridScaleDiff = BasketDiff<Hybrid, FitScaleDer, CovariancePlaneDer>;
    using HybridSpaceDiff = BasketDiff<Hybrid, FitSpaceDer, CovariancePlaneDer>;
    using HybridScaleSpaceDiff = BasketDiff<Hybrid, FitScaleSpaceDer, CovariancePlaneDer>;

    KdTreeDense<Point> tree;
    Scalar scale = generateData(tree);

    for(int i = 0; i < g_repeat; ++i)
    {
        // Single Primitive
        CALL_SUBTEST((testBasicFunctionalities<TestPlane>(tree, scale) ));
        CALL_SUBTEST((testBasicFunctionalities<Sphere>(tree, scale) ));
        // Hybrid
        CALL_SUBTEST((testBasicFunctionalities<Hybrid>(tree, scale) ));
        //  Plane diffs
        CALL_SUBTEST((testBasicFunctionalities<PlaneScaleDiff>(tree, scale) ));
        CALL_SUBTEST((testBasicFunctionalities<PlaneSpaceDiff>(tree, scale) ));
        CALL_SUBTEST((testBasicFunctionalities<PlaneScaleSpaceDiff>(tree, scale) ));
        // Hybrid diffs
        // TODO : Fix hybrid diffs (function calls are ambiguous on primitives)
        // CALL_SUBTEST((testBasicFunctionalities<HybridScaleDiff>(tree, scale) ));
        // CALL_SUBTEST((testBasicFunctionalities<HybridSpaceDiff>(tree, scale) ));
        // CALL_SUBTEST((testBasicFunctionalities<HybridScaleSpaceDiff>(tree, scale) ));

        // Check that we get the same Sphere, whatever the extensions
        auto checkIsSameSphere = [](const auto&f1, const auto&f2){isSameSphere(f1,f2);};
        CALL_SUBTEST((testIsSame<Sphere, Hybrid>(tree, scale, checkIsSameSphere) ));
        CALL_SUBTEST((testIsSame<Sphere, HybridScaleDiff>(tree, scale, checkIsSameSphere) ));
        CALL_SUBTEST((testIsSame<Sphere, HybridSpaceDiff>(tree, scale, checkIsSameSphere) ));
        CALL_SUBTEST((testIsSame<Sphere, HybridScaleSpaceDiff>(tree, scale, checkIsSameSphere) ));

        // Check that we get the same Plane, whatever the extensions
        auto checkIsSamePlane = [](const auto&f1, const auto&f2){isSamePlane(f1,f2);};
        CALL_SUBTEST((testIsSame<TestPlane, Hybrid>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<TestPlane, HybridScaleDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<TestPlane, HybridSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<TestPlane, HybridScaleSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneScaleDiff, HybridSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneScaleDiff, HybridScaleSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneSpaceDiff, HybridScaleDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneSpaceDiff, HybridScaleSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneScaleSpaceDiff, HybridScaleDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneScaleSpaceDiff, HybridSpaceDiff>(tree, scale, checkIsSamePlane) ));
        CALL_SUBTEST((testIsSame<PlaneScaleSpaceDiff, HybridScaleSpaceDiff>(tree, scale, checkIsSamePlane) ));

        auto checkIsSamePlaneDerivative = [](const auto&f1, const auto&f2){
            isSamePlane(f1,f2);
            hasSamePlaneDerivatives(f1, f2);
        };
        // Check that we get the same Plane derivative, whatever if we have one or more primitive fitted.
        int fff = 0;
        CALL_SUBTEST((testIsSame<PlaneScaleDiff, HybridScaleDiff>(tree, scale, checkIsSamePlaneDerivative) ));
        CALL_SUBTEST((testIsSame<PlaneSpaceDiff, HybridSpaceDiff>(tree, scale, checkIsSamePlaneDerivative) ));
        CALL_SUBTEST((testIsSame<PlaneScaleSpaceDiff, HybridScaleSpaceDiff>(tree, scale, checkIsSamePlaneDerivative) ));
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
