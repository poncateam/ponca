/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cuda/ponca_fit_kdtree.cu
 * \brief Example that uses the Fitting and SpatialPartitioning module with Sycl
 * \authors Bastien DOIGNIES
 */

#include <Ponca/Ponca>
#include <sycl/sycl.hpp>
#include "sycl_kdtree.h"

template <typename Fit, typename DataPoint, typename GPUBuffers>
typename DataPoint::Scalar test_fit(GPUBuffers buffers, int i, typename DataPoint::Scalar scale)
{
    KdTreeGPU<DataPoint> kdtree(buffers);

    typename DataPoint::VectorType pos = kdtree.points()[i].pos();

    Fit fit;
    fit.setNeighborFilter({pos, scale});
    fit.computeWithIds(kdtree.rangeNeighbors(i, scale), kdtree.points());

    if (!fit.isStable())
        return NAN;
    return fit.potential(pos);
}

template <typename Scalar, unsigned int Dim>
void testPlaneSycl(sycl::device device, const bool _bUnoriented = false, const bool _bAddPositionNoise = false,
                   const bool _bAddNormalNoise = false)
{
    std::cout << "Running on: " << device.get_info<sycl::info::device::name>() << std::endl;

    using DataPoint        = Ponca::PointPositionNormal<Scalar, Dim>;
    using WeightSmoothFunc = Ponca::DistWeightFunc<DataPoint, Ponca::SmoothWeightKernel<Scalar>>;
    using MeanFitSmooth    = Ponca::Basket<DataPoint, WeightSmoothFunc, Ponca::MeanPlaneFit>;
    using VectorType       = typename DataPoint::VectorType;

    // Point cloud parameters for the plane
    const unsigned int nbPoints = Eigen::internal::random<int>(100, 1000);
    const Scalar width          = Eigen::internal::random<Scalar>(1., 10.);
    const Scalar height         = width;
    const Scalar analysisScale  = Scalar(15.) * std::sqrt(width * height / nbPoints);
    const Scalar centerScale    = Eigen::internal::random<Scalar>(1, 10000);
    const VectorType center     = VectorType::Random() * centerScale;
    const VectorType direction  = VectorType::Random().normalized();

    // Generate the point cloud
    std::vector<DataPoint> points(nbPoints);
    std::vector<Scalar> potentials(nbPoints);
    for (unsigned int i = 0; i < nbPoints; ++i)
    {
        potentials[i] = NAN;
        points[i] = Ponca::getPointOnPlane<DataPoint>(center, direction, width, _bAddPositionNoise, _bAddNormalNoise,
                                                      _bUnoriented);
    }

    //! [Build KdTree on CPU]
    Ponca::KdTreeDense<DataPoint> kdtree(points);
    //! [Build KdTree on CPU]

    auto queue = sycl::queue(device);
    // The size of the data we send between Host and Device
    const unsigned long scalarBufferSize = nbPoints * sizeof(Scalar);
    const unsigned long vectorBufferSize = scalarBufferSize * Dim;

    //! [Copy KdTree on GPU]
    using BuffersGPU = typename KdTreeGPU<DataPoint>::Buffers;

    BuffersGPU deviceKDTree; // Host Buffers referencing data on the device, used to free memory
    deepCopyBuffersToDevice<Ponca::KdTreePointerTraits<DataPoint>>(queue, kdtree.buffers(), deviceKDTree);
    //! [Copy KdTree on GPU]

    // Allocate and copy data on host
    DataPoint* devicePoints  = sycl::malloc_device<DataPoint>(nbPoints, queue);
    Scalar* devicePotentials = sycl::malloc_device<Scalar>(nbPoints, queue);

    queue.memcpy(devicePoints, points.data(), nbPoints * sizeof(DataPoint));
    queue.memcpy(devicePotentials, potentials.data(), nbPoints * sizeof(Scalar));
    queue.wait(); // Wait for "everything", which is just copies

    auto env = queue.submit([&](sycl::handler& h) {
        h.parallel_for(sycl::range<1>(nbPoints), [=](sycl::id<1> i) {
            size_t idx            = i[0];
            devicePotentials[idx] = test_fit<MeanFitSmooth, DataPoint>(deviceKDTree, idx, analysisScale);
        });
    });

    // Get back data
    queue.memcpy(potentials.data(), devicePotentials, nbPoints * sizeof(Scalar), env).wait();

    // Free data
    sycl::free(devicePotentials, queue);
    sycl::free(devicePoints, queue);
}

int main(const int /*argc*/, char** /*argv*/)
{
    // For AdaptiveCpp: we can control what device the code will run with the
    // environment variable ACPP_DEFAULT_SELECTOR_BEHAVIOR.
    testPlaneSycl<double, 3>(sycl::device(sycl::default_selector_v));
    return 0;
}
