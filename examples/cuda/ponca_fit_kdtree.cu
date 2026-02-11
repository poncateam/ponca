#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/Common/pointTypes.h>
#include <Ponca/src/Common/pointGeneration.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTreeTraits.h>
#include <iostream>

#include "cuda_utils.cu"

/*! \brief Computes the fitting process for each point of the point cloud and returns the potential and primitive gradient result.
 *
 * \tparam DataPoint The DataPoint type.
 * \tparam Fit The Fit that will be computed by the Kernel.
 * \param points As an Input, the point cloud.
 * \param nbPoints The number of points in the point cloud
 * \param indices The indices buffer
 * \param nbIndices The number of indices
 * \param nodes The KdTree nodes (internal to the kdtree)
 * \param nbNodes The number of nodes
 * \param analysisScale The radius of the neighborhood.
 * \param potentialResults As an Output, the potential results of the fit for each point of the Point Cloud.
 * \param gradientResults As an Output, the primitiveGradient results of the fit for each point of the Point Cloud.
 */
template<typename DataPoint, typename Fit, typename Node, typename Scalar>
__global__ void fitPotentialAndGradientKernel(
    DataPoint* points,
    const size_t nbPoints,
    int* indices,
    const size_t nbIndices,
    Node * nodes,
    const size_t nbNodes,
    const Scalar analysisScale,
    Scalar* const potentialResults,
    Scalar* const gradientResults
) {
    using VectorType = typename DataPoint::VectorType;

    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

    // Skip when not in the point cloud
    if (i >= nbPoints) return;

    VectorType pos = points[i].pos();
    KdTreeGPU<DataPoint> kdtree(points, nodes, indices, nbPoints, nbNodes, nbIndices);

    // Set up the fit
    Fit fit;
    fit.setNeighborFilter({ pos, analysisScale });

    // Computes the fit
    fit.init();
    for (int j : kdtree.rangeNeighbors(i, analysisScale)) {
        fit.addNeighbor( points[j] );
    }
    fit.finalize();

    // Returns NaN if not stable
    if (! fit.isStable()) {
        potentialResults[i] = NAN;
        for (int d = 0; d < DataPoint::Dim; ++d) {
            gradientResults[i*DataPoint::Dim + d] = NAN;
        }
        return;
    }

    // Return the fit.potential result as an output
    potentialResults[i]   = fit.potential(pos);
    const VectorType grad = fit.primitiveGradient(pos);
    for (int d = 0; d < DataPoint::Dim; ++d) {
        gradientResults[i*DataPoint::Dim + d] = grad(d);
    }
}

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
    typedef Ponca::PointPositionNormalLateBinding<Scalar, Dim> DataPoint;
    typedef Ponca::DistWeightFunc<DataPoint, Ponca::SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef Ponca::Basket<DataPoint, WeightSmoothFunc, Ponca::MeanPlaneFit> MeanFitSmooth;
    typedef typename DataPoint::VectorType VectorType;

    // Point cloud parameters for the plane
    const int nbPoints  = Eigen::internal::random<int>(100, 1000);
    const Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    const Scalar height = width;
    const Scalar analysisScale = Scalar(15.) * std::sqrt( width * height / nbPoints);
    const Scalar centerScale   = Eigen::internal::random<Scalar>(1, 10000);
    const VectorType center    = VectorType::Random() * centerScale;
    const VectorType direction = VectorType::Random().normalized();

    // Interlaced array of position and normal values
    auto* const interlacedArray = new Scalar[nbPoints * Dim * 2];

    for(unsigned int i = 0; i < nbPoints; ++i) {
        auto point = Ponca::getPointOnPlane<Ponca::PointPositionNormal<Scalar, Dim>>(
            center, direction, width,
            _bAddPositionNoise, _bAddNormalNoise, _bUnoriented
        );
        // Fill the interlaced array with the positions and normals.
        for (int d = 0; d<Dim; ++d) {
            interlacedArray[i * Dim * 2 + d]       = point.pos()(d);
            interlacedArray[i * Dim * 2 + d + Dim] = point.normal()(d);
        }
    }

    // Make a point buffer linking to the interlaced array
    std::vector<DataPoint> points;
    points.reserve(nbPoints);
    for (int i = 0; i < nbPoints; ++i) {
        points.emplace_back(interlacedArray, i);
    }

    // Send the internal buffers of the KdTree to the GPU
    Ponca::KdTreeDense<DataPoint> kdtree(points);

    typedef typename Ponca::KdTreeDefaultTraits<DataPoint>::NodeType NodeType;

    std::cout << "nodes size : " << kdtree.nodeCount() << std::endl;

    // The size of the data we send between Host and Device
    const unsigned long pointBufferSize      = nbPoints * sizeof(DataPoint);
    const unsigned long scalarBufferSize     = nbPoints * sizeof(Scalar);
    const unsigned long vectorBufferSize     = scalarBufferSize * Dim;
    const unsigned long interlacedBufferSize = vectorBufferSize * 2;
    const unsigned long nodeBufferSize       = kdtree.nodeCount() * sizeof(NodeType);
    const unsigned long indexBufferSize      = kdtree.sampleCount() * sizeof(int);


    // Send inputs to the GPU (Host to Device)
    Scalar* interlacedArrayDevice;
    DataPoint* pointsDevice;
    int* indicesDevice;
    NodeType* nodesDevice;

    CUDA_CHECK(cudaMalloc(&interlacedArrayDevice, interlacedBufferSize ));
    CUDA_CHECK(cudaMalloc(&pointsDevice         , pointBufferSize ));
    CUDA_CHECK(cudaMalloc(&indicesDevice        , indexBufferSize));
    CUDA_CHECK(cudaMalloc(&nodesDevice          , nodeBufferSize));

    CUDA_CHECK(cudaMemcpy(interlacedArrayDevice , interlacedArray            , interlacedBufferSize, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(pointsDevice          , points.data()          , pointBufferSize     , cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(indicesDevice         , kdtree.samples().data(), indexBufferSize     , cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(nodesDevice           , kdtree.nodes().data()  , nodeBufferSize      , cudaMemcpyHostToDevice));

    // Prepare output buffers
    auto* const potentialResults = new Scalar[nbPoints];
    auto* const gradientResults  = new Scalar[nbPoints*Dim];
    Scalar* potentialResultsDevice;
    Scalar* gradientResultsDevice;
    CUDA_CHECK(cudaMalloc(&potentialResultsDevice, scalarBufferSize));
    CUDA_CHECK(cudaMalloc(&gradientResultsDevice , vectorBufferSize));

    // Set block and grid size depending on number of points
    constexpr int blockSize = 128;
    const     int gridSize  = (nbPoints + blockSize - 1) / blockSize;

    // Binds the points to their internal values on the device
    bindPointsKernel<DataPoint><<<gridSize, blockSize>>>(pointsDevice, interlacedArrayDevice, nbPoints);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Compute the fitting in the kernel
    fitPotentialAndGradientKernel<DataPoint, MeanFitSmooth><<<gridSize, blockSize>>>(
        pointsDevice , nbPoints,
        indicesDevice, kdtree.sampleCount(),
        nodesDevice  , kdtree.nodeCount(),
        analysisScale,
        potentialResultsDevice, gradientResultsDevice
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Fetch the results (Device to Host)
    CUDA_CHECK(cudaMemcpy(potentialResults, potentialResultsDevice, scalarBufferSize, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gradientResults , gradientResultsDevice , vectorBufferSize, cudaMemcpyDeviceToHost));

    // Free CUDA memory
    CUDA_CHECK(cudaFree(pointsDevice));
    CUDA_CHECK(cudaFree(interlacedArrayDevice));
    CUDA_CHECK(cudaFree(potentialResultsDevice));
    CUDA_CHECK(cudaFree(gradientResultsDevice));

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

    delete[] interlacedArray;
    delete[] potentialResults;
    delete[] gradientResults;
}


__host__ int main(const int /*argc*/, char** /*argv*/) {
    std::cout << "Test plane fitting on CUDA..." << std::endl;
    testPlaneCuda<float, 3>();
    std::cout << "(ok)" << std::endl;
}
