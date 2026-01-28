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
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>

#include <vector>

using namespace std;
using namespace Ponca;

template<typename DataPoint>
typename DataPoint::Scalar generateData(vector<DataPoint>& points)
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

    points.reserve(nbPoints);
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < nbPoints; ++i) {
        points.push_back(getPointOnSphere<DataPoint>(radius, center, false, false, false));
    }

    return analysisScale;
}

/*! \brief Variadic function to copy a vector of points to one or multiple vectors of points that are bound to a single interlaced array of positions and normals.
 *
 * \tparam DataPoint Regular point data type (e.g. \ref PointPositionNormal)
 * \tparam DataPointRef Point data types referencing the interlaced array (e.g. \ref PointPositionNormalBinding, \ref PointPositionNormalLateBinding)
 * \param points As an input, a vector of DataPoint.
 * \param pointsBinding As an output, empty vectors of DataPointRef types.
 * \return An interlaced array containing the position and normal values. The `pointsBinding` vectors are linked to this array.
 */
template<typename DataPoint, typename... DataPointRef>
typename DataPoint::Scalar * copyDataRef(std::vector<DataPoint>& points, std::vector<DataPointRef>&... pointsBinding)
{
    using VectorType            = typename DataPoint::VectorType;
    using Scalar                = typename DataPoint::Scalar;

    constexpr int DIM           =  DataPoint::Dim;
    const int nPoints           = int(points.size());
    auto* const interlacedArray = new Scalar[2*DIM*nPoints];
    (pointsBinding.reserve(nPoints), ...);

    for(int i=0; i<nPoints; ++i)
    {
        // We use Eigen Vectors to compute both coordinates and normals,
        // and then copy the raw values to an interlaced array.
        VectorType n = points[i].normal();
        VectorType p = points[i].pos();

        // Grab coordinates and store them as raw buffer
        memcpy(interlacedArray+2*DIM*i    , p.data(), DIM*sizeof(Scalar));
        memcpy(interlacedArray+2*DIM*i+DIM, n.data(), DIM*sizeof(Scalar));

        (pointsBinding.emplace_back(interlacedArray, i), ...);
    }

    return interlacedArray;
}

/*! \brief Verify that the fit results between two point data types are the same. The points are passed through an spatial partitioning structure
 *
 * The fit is computed for each point of the point cloud. The results are compared using the `fit.isApprox` method
 *
 * \tparam Fit A fit class that performs the fitting over the neighborhood
 * \tparam SpatialStruct1 The first Spatial partitioning structure containing a rangeNeighbors method (e.g. \ref KdTree or \ref KnnGraph)
 * \tparam SpatialStruct2 The second Spatial partitioning structure containing a rangeNeighbors method (e.g. \ref KdTree or \ref KnnGraph)
 * \param spatialStruct1 A spatial partitioning structure containing the first set of data point
 * \param spatialStruct2 A spatial partitioning structure the first set of data point that we are comparing to the first
 * \param analysisScale The radius of the neighborhood for the fitting process
 */
template<template<typename> typename Fit, typename SpatialStruct1, typename SpatialStruct2>
void compareFitOverPointTypes( SpatialStruct1& spatialStruct1, SpatialStruct2& spatialStruct2, typename SpatialStruct2::DataPoint::Scalar analysisScale)
{
    using DataPoint1          = typename SpatialStruct1::DataPoint;
    using DataPoint2          = typename SpatialStruct2::DataPoint;
    constexpr int  DIMENSION  = DataPoint1::Dim;
    static_assert( DIMENSION == DataPoint2::Dim, "Both dimension should be the same" );
    static_assert( std::is_same_v<typename DataPoint1::Scalar, typename DataPoint2::Scalar>, "Both scalar type should be the same" );
    const std::vector<DataPoint1>& points1 = spatialStruct1.points();
    const std::vector<DataPoint2>& points2 = spatialStruct2.points();
    VERIFY( points1.size() == points2.size() );

    // Quick testing is requested for coverage
    const int nPoint = QUICK_TESTS ? 1 : int(points1.size());
    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < nPoint; ++i)
    {
        Fit<DataPoint1> f1;
        f1.setNeighborFilter({points1[i].pos(), analysisScale});
        const auto neighborhoodRange1 = spatialStruct1.rangeNeighbors(points1[i].pos(), analysisScale);
        f1.computeWithIds( neighborhoodRange1, points1 );

        Fit<DataPoint2> f2;
        f2.setNeighborFilter({points2[i].pos(), analysisScale});
        const auto neighborhoodRange2 = spatialStruct2.rangeNeighbors(points2[i].pos(), analysisScale);
        f2.computeWithIds( neighborhoodRange2, points2 );

        if (!f1.isStable() || !f2.isStable())
            continue;

        VERIFY(f1.isApprox(f2));
        VERIFY(f2.isApprox(f1));
    }
}

//! \brief Smooth weight neighbor filter class templated over the point data type.
template <typename DataPoint>
using NeighborFilter = DistWeightFunc<DataPoint, SmoothWeightKernel<typename DataPoint::Scalar> >;
//! \brief Fitting templated over the point data type.
template <typename DataPoint>
using TestSphereFit  = Basket<DataPoint, NeighborFilter<DataPoint>, OrientedSphereFit>;

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;
    typedef PointPositionNormalBinding<Scalar, Dim> PointRef;
    typedef PointPositionNormalLateBinding<Scalar, Dim> PointLateRef;

    for(int i = 0; i < g_repeat; ++i)
    {
        // Points to compare
        vector<Point>        points;
        vector<PointRef>     pointsRef;
        vector<PointLateRef> pointsLateRef;

        const Scalar  analysisScale   = generateData(points);
        // Copy the point data to an external buffer, and bind some point vectors to it.
        const Scalar* interlacedArray = copyDataRef(points, pointsRef, pointsLateRef);

        KdTreeDense<Point>        kdtree(points);
        KdTreeDense<PointRef>     kdtreeRef(pointsRef);
        KdTreeDense<PointLateRef> kdtreeLateRef(pointsLateRef);

        // Compare fits made with the kdtree
        CALL_SUBTEST((compareFitOverPointTypes<TestSphereFit>(kdtree, kdtreeRef    , analysisScale)));
        CALL_SUBTEST((compareFitOverPointTypes<TestSphereFit>(kdtree, kdtreeLateRef, analysisScale)));

        // Delete buffer before the next pass
        delete[] interlacedArray;
    }
}

int main(const int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Binding point type in 3 dimensions: float" << flush;
    CALL_SUBTEST_1((callSubTests<float, 3>()));
    cout << " (ok), double" << flush;
    CALL_SUBTEST_2((callSubTests<double, 3>()));
    cout << " (ok)" << flush;
    cout << ", long double" << flush;
    CALL_SUBTEST_3((callSubTests<long double, 3>()));
    cout << " (ok)" << flush;
}
