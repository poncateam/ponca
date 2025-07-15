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

template<typename Fit>
void testBasicFunctionalities(const KdTree<typename Fit::DataPoint>& tree, typename Fit::Scalar analysisScale)
{

    using DataPoint = typename Fit::DataPoint;

    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename Fit::WFunctor WeightFunc;

    const auto& vectorPoints = tree.points();
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
    using fit_CNC = CNC<Point, NoWeightFunc<Point>>;
    //! [CNCFitType]

    // KdTreeDense<Point> tree;
    // Scalar scale = generateData(tree);
    constexpr int N = 100;
	using VectorContainer = typename KdTreeSparse<Point>::PointContainer;
    auto points = VectorContainer(N);
    typedef typename Point::VectorType VectorType;

    std::generate(points.begin(), points.end(), []() {return Point(VectorType::Random()); });

    /// [KdTree pointer usage]
    // Abstract pointer type that can receive KdTreeSparse or KdTreeDense objects
    KdTree<Point> *kdtree {nullptr};
    /// [KdTree pointer usage]


    std::vector<int> sampling;
    sampling.resize(N);
    std::iota(sampling.begin(), sampling.end(), 0);
    /// [KdTree assign dense]
    kdtree = new KdTreeDense<Point> (points);
    /// [KdTree assign dense]

    // CALL_SUBTEST((testBasicFunctionalities<TestPlane>(kdtree, 1) ));
}

int main(const int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test for the CorrectedNormalCurrent fit method" << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
