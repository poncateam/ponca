/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/src/cnc.cpp
 \brief Test validity of the CNC curvature estimator procedure
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/cnc.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>

#include <vector>

#include "Ponca/src/Fitting/curvature.h"
#include "Ponca/src/Fitting/curvatureEstimation.h"
#include "Ponca/src/Fitting/mlsSphereFitDer.h"
#include "Ponca/src/Fitting/weightKernel.h"


using namespace std;
using namespace Ponca;

template<typename DataPoint>
typename DataPoint::Scalar generateSpherePC(KdTree<DataPoint>& tree) {
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
    for(int i = 0; i < int(vectorPoints.size()); ++i) {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, false, false, false);
    }

    tree.clear();
    tree.build(vectorPoints);

    return analysisScale;
}

template<typename Fit>
void testBasicFunctionalities(const KdTree<typename Fit::DataPoint>& tree) {
    const auto& vectorPoints = tree.points();
    auto rng = std::default_random_engine {};
    typename Fit::Scalar eps = testEpsilon<typename Fit::Scalar>()*2;
    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for (int i = 0; i < static_cast<int>(vectorPoints.size()); ++i) {
        const auto &fitInitPoints = vectorPoints[i];

        //! [Fit compute]
        Fit fit1;
        fit1.setEvalPoint(fitInitPoints);
        fit1.compute( tree.points() );
        //! [Fit compute]
        VERIFY(fit1 == fit1);
        VERIFY(! (fit1 != fit1));

        // compute the indices list
        std::vector<int> pointsIndex;
        for (int j = 0; j < tree.points().size(); j++) {
            pointsIndex.push_back(j);
        }
        // Shuffling the indices shouldn't change the outcome of this test
        std::shuffle(std::begin(pointsIndex), std::end(pointsIndex), rng);

        //! [Fit computeWithIds]
        Fit fit2;
        fit2.setEvalPoint(fitInitPoints);
        fit2.computeWithIds( pointsIndex, tree.points() );
        //! [Fit computeWithIds]

        VERIFY((fit2 == fit2));
        VERIFY(! (fit2 != fit2));

        // Compare computeWithIds with compute result
        VERIFY((fit1.isApprox(fit2, eps)));
        VERIFY((fit2.isApprox(fit1, eps)));
    }
}

/// \breif Compare the GaussianCurvature and kMean between two fit
template<typename Fit1, typename Fit2>
void testCompareFit(const KdTree<typename Fit1::DataPoint>& tree, typename Fit1::Scalar analysisScale) {
    typedef typename Fit1::DataPoint DataPoint;
    typedef typename DataPoint::Scalar Scalar;
    const auto& vectorPoints = tree.points();

    // Test for each point if the curvature results are equivalent
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i) {
        // compute the indices list
        std::vector<int> pointsIndex;
        for (const int j : tree.range_neighbors(i, analysisScale)) {
            pointsIndex.push_back(j);
        }

        Fit1 fit1;
        fit1.setWeightFunc({vectorPoints[i].pos(), analysisScale});
        fit1.computeWithIds(pointsIndex, vectorPoints);

        Fit2 fit2;
        fit2.setEvalPoint(vectorPoints[i]);
        fit2.computeWithIds(pointsIndex, vectorPoints);

        Scalar eps = testEpsilon<Scalar>();

        // Compare Fit1 with Fit2
        VERIFY(((fit1.kMean() - fit2.kMean()) < eps));
        VERIFY(((fit1.GaussianCurvature() - fit2.GaussianCurvature()) < eps));
    }
}

template<typename Scalar, int Dim>
void callSubTests() {
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef typename Point::VectorType VectorType;
    //! [SpecializedPointType]
    using SmoothWeightFunc   = Ponca::DistWeightFunc<Point, SmoothWeightKernel<Scalar> >;
    using FitASODiff = BasketDiff<
            Basket<Point, SmoothWeightFunc, Ponca::OrientedSphereFit>,
            DiffType::FitSpaceDer,
            OrientedSphereDer, Ponca::MlsSphereFitDer,
            CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    //! [CNCFitType]
    using Fit_CNC_Independent = CNC<Point, TriangleGenerationMethod::IndependentGeneration>;
    using Fit_CNC_Uniform = CNC<Point, TriangleGenerationMethod::UniformGeneration>;
    using Fit_CNC_Hexagram = CNC<Point, TriangleGenerationMethod::HexagramGeneration>;
    using Fit_CNC_AvgHexagram = CNC<Point, TriangleGenerationMethod::AvgHexagramGeneration>;
    //! [CNCFitType]

    KdTreeDense<Point> tree;
    Scalar analysisScale = generateSpherePC(tree);
    CALL_SUBTEST((testBasicFunctionalities<Fit_CNC_Independent>(tree) ));
    CALL_SUBTEST((testBasicFunctionalities<Fit_CNC_Uniform>(tree) ));
    CALL_SUBTEST((testBasicFunctionalities<Fit_CNC_Hexagram>(tree) ));

    // Compare with ASO
    CALL_SUBTEST((testCompareFit<FitASODiff, Fit_CNC_Independent>(tree, analysisScale) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, Fit_CNC_Uniform>(tree, analysisScale) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, Fit_CNC_Hexagram>(tree, analysisScale) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, Fit_CNC_AvgHexagram>(tree, analysisScale) ));
}

int main(const int argc, char** argv) {
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test for the CorrectedNormalCurrent fit method" << endl;

    cout << "Tests CNC functions in 3 dimensions: float" << flush;
    callSubTests<float, 3>();
    cout << " (ok), double" << flush;
    callSubTests<double, 3>();
    cout << " (ok), long double" << flush;
    callSubTests<long double, 3>();
    cout << " (ok)" << flush;
}
