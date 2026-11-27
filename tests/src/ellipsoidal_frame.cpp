/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/ellipsoidal_frame.cpp
 * \brief Test validity of EllipsoidalNeighborhoodFrame
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/Fitting>

#include <utility>
#include <vector>
#include <cmath>

using namespace std;

template <typename DataPoint, typename WeightKernel>
void testDistWeightFilter()
{
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using MatrixType     = typename DataPoint::MatrixType;
    using Filter         = Ponca::DistWeightEllipsoidalFilter<DataPoint, WeightKernel>;
    using CenteredFilter = Ponca::DistWeightFilter<DataPoint, WeightKernel>;
    static_assert(Ponca::IsNeighborhoodFilter<Filter>);
    static_assert(Ponca::IsNeighborhoodFrame<Filter>);

    const Scalar epsilon     = testEpsilon<Scalar>();
    const Scalar centerScale = Eigen::internal::random<Scalar>(1, 5);
    const VectorType center  = centerScale * VectorType::Random();
    const Scalar t           = Eigen::internal::random<Scalar>(Scalar(0.5), Scalar(3));

    // Verify identity matrix corresponds to centered frame
    {
        const Filter identityFilter(center, t);
        const CenteredFilter centeredFilter(center, t);

        const VectorType q = identityFilter.convertToGlobalBasis(VectorType::Random().normalized() * t);
        const DataPoint dq(q);

        // operator()
        VERIFY(std::abs(identityFilter(dq).first - centeredFilter(dq).first) < epsilon);
        VERIFY((identityFilter(dq).second - centeredFilter(dq).second).norm() < epsilon);

        // Space derivatives
        VERIFY((identityFilter.spacedw(q, dq) - centeredFilter.spacedw(q, dq)).norm() < epsilon);
        VERIFY((identityFilter.spaced2w(q, dq) - centeredFilter.spaced2w(q, dq)).norm() < epsilon);
        // Scale derivatives
        VERIFY(std::abs(identityFilter.scaledw(q, dq) - centeredFilter.scaledw(q, dq)) < epsilon);
        VERIFY(std::abs(identityFilter.scaled2w(q, dq) - centeredFilter.scaled2w(q, dq)) < epsilon);
        // Cross derivatives
        VERIFY((identityFilter.scaleSpaced2w(q, dq) - centeredFilter.scaleSpaced2w(q, dq)).norm() < epsilon);
    }

    // Verify that the distances for points on a sphere with centered filter is the same as point on a ellipse
    // with an ellipsoidal filter
    {
        // Rotation: Q part of QR decomposition
        const MatrixType rotation = MatrixType::Random().householderQr().householderQ();
        VectorType scales;
        for (int i = 0; i < DataPoint::Dim; ++i)
            scales(i) = Eigen::internal::random<Scalar>(Scalar(0.5), Scalar(2));
        const MatrixType axes = rotation * scales.asDiagonal();

        const Filter ellipsoidFilter(center, t, axes);
        const CenteredFilter sphereFilter(center, t);

        const VectorType v          = VectorType::Random();
        const VectorType pSphere    = center + v;
        const VectorType pEllipsoid = center + axes * v;
        const DataPoint dSphere(pSphere);
        const DataPoint dEllipsoid(pEllipsoid);

        // Weights
        VERIFY(std::abs(ellipsoidFilter(dEllipsoid).first - sphereFilter(dSphere).first) < epsilon);

        // Scale derivatives should also match as it depends only on distances, space depends on position
        VERIFY(std::abs(ellipsoidFilter.scaledw(pEllipsoid, dEllipsoid) - sphereFilter.scaledw(pSphere, dSphere)) <
               epsilon);
        VERIFY(std::abs(ellipsoidFilter.scaled2w(pEllipsoid, dEllipsoid) - sphereFilter.scaled2w(pSphere, dSphere)) <
               epsilon);
    }
}

template <typename Scalar, int Dim>
void callSubTests()
{
    using Point = Ponca::PointPositionNormal<Scalar, Dim>;

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST((testDistWeightFilter<Point, Ponca::SmoothWeightKernel<Scalar>>()));
        CALL_SUBTEST((testDistWeightFilter<Point, Ponca::GaussianWeightKernel<Scalar>>()));
    }
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test EllipsoidalNeighborhoodFrame : " << flush;

    cout << "float" << flush;
    callSubTests<float, 3>();
    cout << " (ok), " << flush;
    cout << "double" << flush;
    callSubTests<double, 3>();
    callSubTests<double, 2>();
    cout << " (ok), " << flush;
    cout << "long double" << flush;
    callSubTests<long double, 3>();
    cout << " (ok)" << endl;
}
