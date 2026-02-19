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

template<typename DataPoint, typename VectorType>
typename DataPoint::Scalar generateSpherePC(
    KdTree<DataPoint>& tree,
    const int nbPoints = Eigen::internal::random<int>(500, 1000),
    const VectorType& center = VectorType::Random() * Eigen::internal::random<typename DataPoint::Scalar>(1,10000)
) {
    typedef typename DataPoint::Scalar Scalar;

    Scalar radius = Eigen::internal::random<Scalar>(1., 10.);
    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
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

/*!
 * \breif Test the immutability of the compute methods
 *
 * Checks if the order of the data matters for the outputs or not by comparing two compute :
 * One with a compute shuffled indices list
 *
 * \tparam Fit Fit type that will be used to compute the elements
 * \tparam Scalar Scalar type
 * \param tree The KdTree
 * \param analysisScale The size of the neighborhood
 * \param epsilon The precision of the mutability check
 */
template<typename Fit, typename Scalar>
void testBasicFunctionalities(
    const KdTree<typename Fit::DataPoint>& tree,
    const Scalar analysisScale,
    const Scalar epsilon = testEpsilon<Scalar>() *2
) {
    using DataPoint = typename Fit::DataPoint;

    const auto& vectorPoints = tree.points();
    auto rng = std::default_random_engine {};

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for (int i = 0; i < static_cast<int>(vectorPoints.size()); ++i) {
        const DataPoint& evalPoint = vectorPoints[i];

        //! [Fit compute]
        Fit fit1;
        fit1.setNeighborFilter({evalPoint, analysisScale});
        fit1.compute( vectorPoints );
        //! [Fit compute]
        VERIFY(fit1 == fit1);
        VERIFY(! (fit1 != fit1));

        // Sample the neighbors
        std::vector<int> pointsIndex;
        pointsIndex.push_back(i);
        for (int j : tree.rangeNeighbors(i, analysisScale)) {
            pointsIndex.push_back(j);
        }

        // Shuffles the indices shouldn't change the outcome of this test
        std::shuffle(std::begin(pointsIndex), std::end(pointsIndex), rng);

        //! [Fit computeWithIds]
        Fit fit2;
        fit2.setNeighborFilter({evalPoint, analysisScale});
        fit2.computeWithIds( pointsIndex, vectorPoints );
        //! [Fit computeWithIds]

        VERIFY((fit1.isStable()));
        VERIFY((fit2.isStable()));

        // Equal to self test
        VERIFY((fit2 == fit2));
        VERIFY(! (fit2 != fit2));
        VERIFY((fit1 == fit1));
        VERIFY(! (fit1 != fit1));

        if (std::abs(fit1.kMean()  - fit2.kMean()) > epsilon)
            std::cout << "std::abs(kMean()  - other.kMean() :" << std::abs(fit1.kMean()  - fit2.kMean()) << std::endl;
        if (std::abs(fit1.GaussianCurvature() - fit2.GaussianCurvature()) > epsilon)
            std::cout << "std::abs(GaussianCurvature() - other.GaussianCurvature() :" << std::abs(fit1.GaussianCurvature() - fit2.GaussianCurvature()) << std::endl;

        // Compare computeWithIds with compute result
        VERIFY((std::abs(fit1.kMean() - fit2.kMean()) < epsilon));
        VERIFY((std::abs(fit1.GaussianCurvature() - fit2.GaussianCurvature()) < epsilon));
    }
}

/// \breif Compare the GaussianCurvature and kMean between two fit
template<typename Fit1, typename Fit2, bool orderedByDistance = false, typename Scalar>
void testCompareFit(
    const KdTree<typename Fit1::DataPoint>& tree,
    const Scalar analysisScale,
    const Scalar epsilon = testEpsilon<Scalar>() *2
) {
    const auto& vectorPoints = tree.points();
    const int nPoints = int(vectorPoints.size());
    // Test for each point if the curvature results are equivalent
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < nPoints; ++i) {
        typename Fit1::NeighborFilter w {vectorPoints[i].pos(), analysisScale};
        // compute the indices list
        std::vector<int> pointsIndex;
        pointsIndex.push_back(i);

        if constexpr (orderedByDistance) {
            for (int j : tree.kNearestNeighbors(i, vectorPoints.size())) {
                // Stops when we go past the analysis scale
                if (w(vectorPoints[ j ]).first == Scalar(0.))
                    break;
                pointsIndex.push_back(j);
            }
        } else {
            for (const int j : tree.rangeNeighbors(i, analysisScale)) {
                pointsIndex.push_back(j);
            }
        }

        Fit1 fit1;
        fit1.setNeighborFilter({vectorPoints[i], analysisScale});
        fit1.compute(vectorPoints);

        Fit2 fit2;
        fit2.setNeighborFilter({vectorPoints[i], analysisScale});
        fit2.computeWithIds(pointsIndex, vectorPoints);

        VERIFY((fit1.isStable()));
        VERIFY((fit2.isStable()));

        // Compare Fit1 with Fit2
        VERIFY((std::abs(fit1.kMean() - fit2.kMean()) < epsilon));
        VERIFY((std::abs(fit1.GaussianCurvature() - fit2.GaussianCurvature()) < epsilon));
    }
}

template<typename Scalar, int Dim>
void callSubTests() {
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef typename Point::VectorType VectorType;
    //! [SpecializedPointType]

    using SmoothWeightFunc   = DistWeightFunc< Point, SmoothWeightKernel<Scalar> >;
    using FitASODiff = BasketDiff<
            Basket< Point, SmoothWeightFunc, OrientedSphereFit >,
            FitSpaceDer,
            OrientedSphereDer, MlsSphereFitDer,
            CurvatureEstimatorBase, NormalDerivativesCurvatureEstimator>;

    //! [CNCFitType]
    using FitCNCIndependent = CNC<Point, IndependentGeneration>;
    using FitCNCUniform     = CNC<Point, UniformGeneration>;
    using FitCNCHexagram    = CNC<Point, HexagramGeneration>;
    using FitCNCAvgHexagram = CNC<Point, AvgHexagramGeneration>;
    //! [CNCFitType]

    // Generate sphere point cloud
    KdTreeDense<Point> tree;
    int nbPoints = QUICK_TESTS
        ? Eigen::internal::random<int>(100 , 200)
        : Eigen::internal::random<int>(5000, 7000) ;
    const VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);
    const Scalar analysisScale = generateSpherePC(tree, nbPoints, center);

    const Scalar highEpsilon {Scalar(0.1)};

    // Tests validity of compute despite index shuffle
    CALL_SUBTEST((testBasicFunctionalities<FitCNCIndependent>(tree, analysisScale, highEpsilon) ));
    CALL_SUBTEST((testBasicFunctionalities<FitCNCUniform>(tree, analysisScale, highEpsilon) ));

    // Compare with ASO
    CALL_SUBTEST((testCompareFit<FitASODiff, FitCNCIndependent>(tree, analysisScale) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, FitCNCUniform>(tree, analysisScale) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, FitCNCHexagram, true>(tree, analysisScale, highEpsilon) ));
    CALL_SUBTEST((testCompareFit<FitASODiff, FitCNCAvgHexagram, true>(tree, analysisScale, highEpsilon) ));
}

int main(const int argc, char** argv) {
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test for the CorrectedNormalCurrent fit method" << endl;

    cout << "Tests CNC functions in 3 dimensions:";
    cout << "float"              << flush;
    callSubTests<float, 3>();
    cout << " (ok), double"      << flush;
    callSubTests<double, 3>();
    cout << " (ok), long double" << flush;
    callSubTests<long double, 3>();
    cout << " (ok)"              << endl;
}
