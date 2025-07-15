/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/src/cnc.cpp
 \brief Test validity of the cnc fitting procedures
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
        fit1.init();
        for(auto it = vectorPoints.begin(); it != vectorPoints.end(); ++it)
            fit1.addNeighbor(*it);
        fit1.finalize();
        //! [Fit Manual Traversal]

        // use compute function
        //! [Fit Compute]
        Fit fit2;
        fit2.compute(vectorPoints);
        //! [Fit Compute]

        // also test comparison operators
        VERIFY(fit1 == fit1);
        VERIFY(fit2 == fit2);
        VERIFY(fit1 == fit2);
        VERIFY(! (fit1 != fit1));
        VERIFY(! (fit1 != fit2));
        VERIFY(! (fit2 != fit2));

        // we skip kdtree test for float: using the kdtree changes the order of the neighbors, which in turn changes the
        // rounding error accumulations, and thus the final result
        if (std::is_same<Scalar, float>::value || std::is_same<Scalar, long double>::value)
            continue;
        //! [Fit computeWithIds]
        Fit fit3;
        fit3.computeWithIds( tree.range_neighbors(fitInitPos, analysisScale), vectorPoints );
        //! [Fit computeWithIds]
        VERIFY(fit3 == fit3);
        VERIFY(fit1 == fit3);
        VERIFY(! (fit1 != fit3));
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    //! [SpecializedPointType]
    typedef PointPositionNormal<Scalar, Dim> Point;
    //! [SpecializedPointType]

    // We test only primitive functions and not the fitting procedure
    //! [WeightFunction]
    using WeightFunc = DistWeightFunc<Point, SmoothWeightKernel<Scalar> >;
    //! [WeightFunction]
    using Sphere     = Basket<Point, WeightFunc, OrientedSphereFit>;
    //! [PlaneFitType]
    using TestPlane = Basket<Point, WeightFunc, CovariancePlaneFit>;
    //! [PlaneFitType]
    //! [FitType]
    using Sphere = Basket<Point, WeightFunc, OrientedSphereFit>;
    //! [FitType]
    //! [CNCFitType]
    using Fit_CNC = CNC<Point, NoWeightFunc<Point>>;
    //! [CNCFitType]

    KdTreeDense<Point> tree;
    Scalar scale = generateData(tree);
    CALL_SUBTEST((testBasicFunctionalities<Fit_CNC>(tree, scale) ));

 //    constexpr int N = 100;
 // using VectorContainer = typename KdTreeSparse<Point>::PointContainer;
 //    auto points = VectorContainer(N);
 //    typedef typename Point::VectorType VectorType;
 //
 //    std::generate(points.begin(), points.end(), []() {return Point(VectorType::Random()); });
 //
 //    /// [KdTree pointer usage]
 //    // Abstract pointer type that can receive KdTreeSparse or KdTreeDense objects
 //    KdTree<Point> *kdtree {nullptr};
 //    /// [KdTree pointer usage]
 //
 //
 //    std::vector<int> sampling;
 //    sampling.resize(N);
 //    std::iota(sampling.begin(), sampling.end(), 0);
 //    /// [KdTree assign dense]
 //    kdtree = new KdTreeDense<Point> (points);
 //    /// [KdTree assign dense]

    // CALL_SUBTEST((testBasicFunctionalities<TestPlane>(kdtree, 1) ));
}

int main(const int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test for the CorrectedNormalCurrent fit method" << endl;

    callSubTests<float, 3>();
    // callSubTests<double, 3>();
    // callSubTests<long double, 3>();
}
