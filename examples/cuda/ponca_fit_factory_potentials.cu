/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cuda/ponca_fit_plane.cu
 * \brief Example that uses the Fitting module with CUDA
 * \authors Auberval Florian, Nicolas Mellado
 */

#include <Ponca/Fitting>
#include <Ponca/src/Common/pointTypes.h>
#include <Ponca/src/Common/pointGeneration.h>
#include <iostream>

#include "cuda_utils.cu"

/*! \brief Computes the fitting process for each point of the point cloud and returns the potential and primitive
 * gradient result.
 *
 * \tparam DataPoint The DataPoint type.
 * \tparam Fit The Fit that will be computed by the Kernel.
 * \param points As an Input, the point cloud.
 * \param nbPoints The total number of points in the point cloud.
 * \param analysisScale The radius of the neighborhood.
 * \param potentialResults As an Output, the potential results of the fit for each point of the Point Cloud.
 */
template <typename DataPoint, typename Fit, typename Scalar>
__global__ void fitPotentialAndGradientKernel(DataPoint* points, const unsigned int nbPoints,
                                              const Scalar analysisScale, Scalar* const potentialResults)
{
    using VectorType = typename DataPoint::VectorType;

    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

    // Skip when not in the point cloud
    if (i >= nbPoints)
        return;

    VectorType pos = points[i].pos();

    // Set up the fit
    Fit fit;
    fit.setNeighborFilter({pos, analysisScale});

    // Computes the fit
    fit.compute(points, points + nbPoints);

    // Returns NaN if not stable
    if (!fit.isStable())
    {
        potentialResults[i] = NAN;
        return;
    }

    // Return the fit.potential result as an output
    potentialResults[i] = fit.potential(pos);
}

/*! \brief Test a MeanPlaneFit on a plane using the CUDA kernel
 *
 * \tparam Scalar The scalar type (e.g. double, float or long double...)
 * \tparam Dim The number of dimension that the VectorType will have.
 * \param _bUnoriented Generates an unoriented point cloud.
 * \param _bAddPositionNoise Determines if we add a randomly generated offset to the position.
 * \param _bAddNormalNoise Determines if we add a randomly generated offset to the normal.
 */
template <typename Scalar, int Dim>
__host__ void testPlaneCuda(const bool _bUnoriented = false, const bool _bAddPositionNoise = false,
                            const bool _bAddNormalNoise = false)
{
    using DataPoint        = Ponca::PointPositionNormal<Scalar, Dim>;
    using WeightSmoothFunc = Ponca::DistWeightFilter<DataPoint, Ponca::SmoothWeightKernel<Scalar>>;
    using VectorType       = typename DataPoint::VectorType;
    using Factory          = Ponca::Factory<DataPoint, WeightSmoothFunc, Ponca::FitSpaceDer>;

    // Point cloud parameters for the plane
    const unsigned int nbPoints = Eigen::internal::random<int>(100, 1000);
    const Scalar width          = Eigen::internal::random<Scalar>(1., 10.);
    const Scalar height         = width;
    const Scalar analysisScale  = Scalar(15.) * std::sqrt(width * height / nbPoints);
    const Scalar centerScale    = Eigen::internal::random<Scalar>(1, 10000);
    const VectorType center     = VectorType::Random() * centerScale;
    const VectorType direction  = VectorType::Random().normalized();

    // The size of the data we send between Host and Device
    const unsigned long pointBufferSize  = nbPoints * sizeof(DataPoint);
    const unsigned long scalarBufferSize = nbPoints * sizeof(Scalar);

    // Generate the point cloud
    std::vector<DataPoint> points(nbPoints);
    for (unsigned int i = 0; i < nbPoints; ++i)
    {
        points[i] = Ponca::getPointOnPlane<DataPoint>(center, direction, width, _bAddPositionNoise, _bAddNormalNoise,
                                                      _bUnoriented);
    }

    // Send inputs to the GPU (Host to Device)
    DataPoint* pointsDevice;
    CUDA_CHECK(cudaMalloc(&pointsDevice, pointBufferSize));
    CUDA_CHECK(cudaMemcpy(pointsDevice, points.data(), pointBufferSize, cudaMemcpyHostToDevice));

    // Prepare output buffers
    auto* const potentialResults = new Scalar[nbPoints];
    Scalar* potentialResultsDevice;
    CUDA_CHECK(cudaMalloc(&potentialResultsDevice, scalarBufferSize));

    // Set block and grid size depending on number of points
    constexpr unsigned int blockSize = 128;
    const unsigned int gridSize      = (nbPoints + blockSize - 1) / blockSize;

    // Compute the fitting in the kernel

    auto list = Factory::template Filter<Ponca::ImplicitPrimitiveProvider>();
    list.foreach ([&] __host__(auto& x) {
        using Fit = decltype(x.object);
        fitPotentialAndGradientKernel<DataPoint, Fit>
            <<<gridSize, blockSize>>>(pointsDevice, nbPoints, analysisScale, potentialResultsDevice);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Fetch the results (Device to Host)
        CUDA_CHECK(cudaMemcpy(potentialResults, potentialResultsDevice, scalarBufferSize, cudaMemcpyDeviceToHost));

        // Validate results
        std::cout << "Results for: " << x.name << ":\n";
        for (int i = 0; i < nbPoints; ++i)
        {
            std::cout << "\t" << i << ": , potential:" << potentialResults[i] << "\n";
        }
    });

    // Free CUDA memory
    CUDA_CHECK(cudaFree(pointsDevice));
    CUDA_CHECK(cudaFree(potentialResultsDevice));

    delete[] potentialResults;
}

__host__ int main(const int /*argc*/, char** /*argv*/)
{
    std::cout << "Test plane fitting on CUDA..." << std::endl;
    testPlaneCuda<float, 3>();
    std::cout << "(ok)" << std::endl;
}
