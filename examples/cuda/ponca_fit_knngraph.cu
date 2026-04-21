/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cuda/ponca_fit_kdtree.cu
 * \brief Example that uses the Fitting and SpatialPartitioning module with CUDA
 * \authors Auberval Florian, Nicolas Mellado
 */


#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/Common/pointTypes.h>
#include <Ponca/src/Common/pointGeneration.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTreeTraits.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>
#include <iostream>

#include "cuda_utils.cu"

template<typename DataPoint>
struct KnnGraphRangeFunctor {
    static __device__ inline auto query(
        KnnGraphGPU<DataPoint>& d_knngraph,
        int i, typename DataPoint::Scalar analysisScale
    ) -> Ponca::KnnGraphRangeQuery<Ponca::KnnGraphPointerTraits<DataPoint>> {
        //! [Use KnnGraph.rangeNeighbors on the GPU]
        return d_knngraph.rangeNeighbors(i, analysisScale);
        //! [Use KnnGraph.rangeNeighbors on the GPU]
    }
};

template<typename DataPoint>
struct KnnGraphKNearestFunctor {
    static __device__ inline auto query(
        KnnGraphGPU<DataPoint>& d_knngraph,
        int i, typename DataPoint::Scalar analysisScale
    ) -> Ponca::KnnGraphKNearestQuery<Ponca::KnnGraphPointerTraits<DataPoint>> {
        //! [Use KnnGraph.kNearestNeighbors on the GPU]
        return d_knngraph.kNearestNeighbors(i);
        //! [Use KnnGraph.kNearestNeighbors on the GPU]
    }
};

/*! \brief Test a MeanPlaneFit on a plane using the CUDA kernel
 *
 * \tparam Scalar The scalar type (e.g. double, float or long double...)
 * \tparam Dim The number of dimension that the VectorType will have.
 * \param _bUnoriented Generates an unoriented point cloud.
 * \param _bAddPositionNoise Determines if we add a randomly generated offset to the position.
 * \param _bAddNormalNoise Determines if we add a randomly generated offset to the normal.
 */
template<typename Scalar, int Dim>
__host__ void testPlaneCuda(
    const bool _bUnoriented       = false,
    const bool _bAddPositionNoise = false,
    const bool _bAddNormalNoise   = false
) {
    using DataPoint        = Ponca::PointPositionNormal<Scalar, Dim>;
    using WeightSmoothFunc = Ponca::DistWeightFunc<DataPoint, Ponca::SmoothWeightKernel<Scalar> >;
    using MeanFitSmooth    = Ponca::Basket<DataPoint, WeightSmoothFunc, Ponca::MeanPlaneFit>;
    using VectorType       = typename DataPoint::VectorType;

    // Point cloud parameters for the plane
    const unsigned int nbPoints = Eigen::internal::random<int>(100, 1000);
    const Scalar width          = Eigen::internal::random<Scalar>(1., 10.);
    const Scalar height         = width;
    const Scalar analysisScale  = Scalar(15.) * std::sqrt( width * height / nbPoints);
    const Scalar centerScale    = Eigen::internal::random<Scalar>(1, 10000);
    const VectorType center     = VectorType::Random() * centerScale;
    const VectorType direction  = VectorType::Random().normalized();

    // Generate the point cloud
    std::vector<DataPoint> points(nbPoints);
    for(unsigned int i = 0; i < nbPoints; ++i) {
        points[i] = Ponca::getPointOnPlane<DataPoint>(
            center, direction, width,
            _bAddPositionNoise, _bAddNormalNoise, _bUnoriented
        );
    }

    //! [Build KnnGraph for CPU]
    Ponca::KdTreeDense<DataPoint> kdtree(points);
    Ponca::KnnGraph<DataPoint>    knngraph(kdtree, points.size()-1);
    //! [Build KnnGraph for CPU]

    std::cout << "Number of nodes in the KdTree : " << kdtree.nodeCount() << std::endl;

    // The size of the data we send between Host and Device
    const unsigned long scalarBufferSize     = nbPoints * sizeof(Scalar);
    const unsigned long vectorBufferSize     = scalarBufferSize * Dim;

    //! [Copy KnnGraph on GPU]
    using BuffersGPU = typename KnnGraphGPU<DataPoint>::Buffers;
    BuffersGPU* knnGraphBuffersDevice;
    CUDA_CHECK(cudaMalloc(&knnGraphBuffersDevice, sizeof(BuffersGPU)));
    BuffersGPU hostBuffersHoldingDevicePointers; // Host Buffers referencing data on the device, used to free memory
    deepCopyKnnGraphBuffersToDevice<Ponca::KdTreePointerTraits<DataPoint>>(
        knngraph.buffers(), hostBuffersHoldingDevicePointers, knnGraphBuffersDevice
    );
    //! [Copy KnnGraph on GPU]

    // Prepare output buffers
    auto* const potentialResults = new Scalar[nbPoints];
    auto* const gradientResults  = new Scalar[nbPoints*Dim];
    Scalar* potentialResultsDevice;
    Scalar* gradientResultsDevice;
    CUDA_CHECK(cudaMalloc(&potentialResultsDevice, scalarBufferSize));
    CUDA_CHECK(cudaMalloc(&gradientResultsDevice , vectorBufferSize));

    // Set block and grid size depending on number of points
    constexpr unsigned int blockSize = 128;
    const     unsigned int gridSize  = (nbPoints + blockSize - 1) / blockSize;

    // Compute the fitting in the kernel
    fitPotentialAndGradientKernel<KnnGraphGPU<DataPoint>, MeanFitSmooth, KnnGraphKNearestFunctor<DataPoint>>
    <<<gridSize, blockSize>>>(
        knnGraphBuffersDevice, analysisScale,           // Inputs
        potentialResultsDevice, gradientResultsDevice   // Outputs
    );

    CUDA_CHECK(cudaGetLastError()); // Catch kernel launch errors
    CUDA_CHECK(cudaDeviceSynchronize());

    // Fetch the results (Device to Host)
    CUDA_CHECK(cudaMemcpy(potentialResults, potentialResultsDevice, scalarBufferSize, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gradientResults , gradientResultsDevice , vectorBufferSize, cudaMemcpyDeviceToHost));

    // Free CUDA memory
    CUDA_CHECK(cudaFree(potentialResultsDevice));
    CUDA_CHECK(cudaFree(gradientResultsDevice));
    //! [Free KnnGraph from memory on GPU]
    freeKnnGraphBuffersOnDevice(hostBuffersHoldingDevicePointers);
    CUDA_CHECK(cudaFree(knnGraphBuffersDevice));
    //! [Free KnnGraph from memory on GPU]

    // Validate results
    const auto epsilon = Scalar(0.001);
    for (int i = 0; i < nbPoints; ++i) {
        const VectorType primGrad = Eigen::Map< const VectorType >(gradientResults + Dim*i  );

        std::cout << "i:" << i << ", potential:" << potentialResults[i] << ", ";
        std::cout << "primitiveGradient:" << primGrad.transpose() << " ; " << std::endl;

        if (std::abs(potentialResults[i]) > epsilon || Scalar(1.) - std::abs(primGrad.dot(direction)) > epsilon)
        {
            std::cerr << "Test failed in " << __FILE__ << " (" << __LINE__ << ")" << std::endl;
            abort();
        }
    }

    delete[] potentialResults;
    delete[] gradientResults;
}


__host__ int main(const int /*argc*/, char** /*argv*/) {
    std::cout << "Example plane fitting using KnnGraph on CUDA..." << std::endl;
    testPlaneCuda<float, 3>();
    std::cout << "(ok)" << std::endl;
}
