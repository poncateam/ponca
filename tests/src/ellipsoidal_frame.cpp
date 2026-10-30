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

template <typename DataPoint>
void testFrame()
{
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    using MatrixType = typename DataPoint::MatrixType;
    using Frame      = Ponca::EllipsoidalNeighborhoodFrame<DataPoint>;
    static_assert(Ponca::IsNeighborhoodFrame<Frame>);

    const Scalar epsilon     = testEpsilon<Scalar>();
    const Scalar centerScale = Eigen::internal::random<Scalar>(1, 5);
    const VectorType center  = centerScale * VectorType::Random();

    // Random orthonormal basis (extract Q from QR factorization of random matrix)
    const MatrixType rotation = MatrixType(MatrixType::Random().householderQr().householderQ());
    
    // Scale column for "distorded" ellipses
    VectorType axes[DataPoint::Dim];
    for (int i = 0; i < DataPoint::Dim; ++i)
        axes[i] = rotation.col(i) * Eigen::internal::random<Scalar>(Scalar(0.5), Scalar(5));

    // Verify the change of basis matrix (i.e. projection of axes to unit vectors)
    const Frame frame(center, axes);
    for (int i = 0; i < DataPoint::Dim; ++i)
    {
        VERIFY((frame.convertToLocalBasis(center + axes[i]) - VectorType::Unit(i)).norm() < epsilon);
    }

    // Verify that covertToGlobal(convertToLocal) = identity
    const VectorType q = VectorType::Random() * Eigen::internal::random<Scalar>(0, 10);
    VERIFY((frame.convertToGlobalBasis(frame.convertToLocalBasis(q)) - q).norm() < epsilon);
    VERIFY((frame.convertToGlobalBasis(frame.convertToLocalBasis(q, false), false) - q).norm() < epsilon);

    // Verify identity matrix corresponds to centered frame
    const Frame identityFrame(center); // Default axes is identity
    const Ponca::CenteredNeighborhoodFrame<DataPoint> centeredFrame(center);
    VERIFY((identityFrame.convertToLocalBasis(q) - centeredFrame.convertToLocalBasis(q)).norm() < epsilon);
}

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
}

template <typename Scalar, int Dim>
void callSubTests()
{
    using Point = Ponca::PointPositionNormal<Scalar, Dim>;

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST((testFrame<Point>()));
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
