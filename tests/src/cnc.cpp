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
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/cnc.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>

#include <vector>


using namespace std;
using namespace Ponca;

template<typename DataPoint>
typename DataPoint::Scalar generateData(KdTree<DataPoint>& tree) {
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
void testBasicFunctionalities(const KdTree<typename Fit::DataPoint>& tree, typename Fit::Scalar analysisScale) {
    const auto& vectorPoints = tree.points();

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        const auto &fitInitPos = vectorPoints[i].pos();

        // use compute function
        //! [Fit Compute]
        Fit fit1;
        fit1.compute(vectorPoints);
        //! [Fit Compute]

        // also test comparison operators
        VERIFY(fit1 == fit1);
        VERIFY(! (fit1 != fit1));

        //! [Fit computeWithIds]
        Fit fit2;
        // Sort fit1
        std::vector<int> neighbors3;
        for (int iNeighbor : tree.range_neighbors(fitInitPos, analysisScale))
            neighbors3.push_back(iNeighbor);
        // Compute the neighbors
        fit2.computeWithIds( neighbors3, vectorPoints );
        //! [Fit computeWithIds]

        VERIFY((fit2 == fit2));
        VERIFY(! (fit2 != fit2));

        // CNC is a stochastic algorithm, so we have to use an epsilon
        // TODO : add an isApprox method
        // VERIFY((fit1.isApprox(fit2)));
        // VERIFY((fit2.isApprox(fit1)));
    }
}

template<typename Scalar, int Dim>
void callSubTests() {
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    //! [SpecializedPointType]

    // We test only primitive functions and not the fitting procedure
    //! [WeightFunction]
    using WeightFunc = DistWeightFunc<Point, SmoothWeightKernel<Scalar> >;
    //! [WeightFunction]
    //! [CNCFitType]
    using Fit_CNC = CNC<Point, NoWeightFunc<Point>>;
    //! [CNCFitType]

    KdTreeDense<Point> tree;
    Scalar scale = generateData(tree);
    CALL_SUBTEST((testBasicFunctionalities<Fit_CNC>(tree, scale) ));
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
