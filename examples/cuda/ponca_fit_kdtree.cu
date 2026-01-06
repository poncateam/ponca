
#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/SpatialPartitioning>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <vector>

#include "Ponca/src/SpatialPartitioning/KdTree/kdTreeTraits.h"
#include "cuda_utils.cu"
#include "../../tests/common/testing.h"
#include "../../tests/common/testUtils.h"


/*! \brief Computes a fit for each point of the point cloud and returns the potential result.
 *
 * \tparam DataPoint The DataPoint type.
 * \tparam Fit The Fit that will be computed by the Kernel.
 * \param positions As an Input, the array of positions of the point cloud.
 * \param normals As an Input, the array of normals of the point cloud.
 * \param nodes
 * \param analysisScale The radius of the neighborhood.
 * \param nbPoints The total number of points in the point cloud.
 * \param potentialResults As an Output, the potential results of the fit for each point of the Point Cloud.
 * \param gradiantResults As an Output, the primitiveGradient results of the fit for each point of the Point Cloud.
 */
template<typename DataPoint, typename Fit, typename Node, typename Scalar>
__global__ void fitPotentialKernel(
    const Scalar* const positions,
    const Scalar* const normals,
    Node * nodes,
    const Scalar analysisScale,
    const int nbPoints,
    Scalar* const potentialResults,
    Scalar* const gradiantResults
) {
    using VectorType = typename DataPoint::VectorType;

    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

    // Skip when not in the point cloud
    if (i >= nbPoints) return;

    // Make the evaluation point of the fit
    const auto evalPoint = makeDataPoint<DataPoint>(i, positions, normals);

    // // Build the kdtree
    // kdtree.prebuild()

    // Set up the fit
    Fit fit;
    fit.setNeighborFilter({ evalPoint.pos(), analysisScale });

    // Computes the fit
    fit.init();
    for (int j = 0; j < nbPoints; ++j) {
        fit.addNeighbor( makeDataPoint<DataPoint>(j, positions, normals) );
    }
    fit.finalize();

    // Returns NaN if not stable
    if (! fit.isStable()) {
        potentialResults[i] = NAN;
        for (int d = 0; d < DataPoint::Dim; ++d) {
            gradiantResults[i*DataPoint::Dim + d] = NAN;
        }
        return;
    }

    // Return the fit.potential result as an output
    potentialResults[i] = fit.potential(evalPoint.pos());
    VectorType grad = fit.primitiveGradient(evalPoint.pos());
    for (int d = 0; d < DataPoint::Dim; ++d) {
        gradiantResults[i*DataPoint::Dim + d] = grad(d);
    }
}

/*! \brief Test a MeanPlaneFit on a plane using the CUDA kernel
 *
 * \tparam Scalar The scalar type (e.g. double, float, long double...)
 * \tparam Dim The number of dimension that the VectorType will have.
 * \param _bUnoriented Generates an unoriented point cloud.
 * \param _bAddPositionNoise Determines if we add a randomly generated offset to the posion.
 * \param _bAddNormalNoise Determines if we add a randomly generated offset to the normal.
 */
template<typename Scalar, int Dim>
__host__ void testPlaneCuda( const bool _bUnoriented = false,  const bool _bAddPositionNoise = false,  const bool _bAddNormalNoise = false)
{
    typedef PointPositionNormalGPU DataPoint;
    typedef Ponca::DistWeightFunc<DataPoint, Ponca::SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef Ponca::Basket<DataPoint, WeightSmoothFunc, Ponca::MeanPlaneFit> MeanFitSmooth;
    typedef typename DataPoint::VectorType VectorType;

    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar height = width;

    Scalar analysisScale = Scalar(15.) * std::sqrt( width * height / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1, 10000);
    VectorType center    = {
        Eigen::internal::random<Scalar>(0., centerScale),
        Eigen::internal::random<Scalar>(0., centerScale),
        Eigen::internal::random<Scalar>(0., centerScale)
    };
    Scalar dx = Eigen::internal::random<Scalar>(0., 1.);
    Scalar dy = Eigen::internal::random<Scalar>(0., 1.);
    Scalar dz = Eigen::internal::random<Scalar>(0., 1.);
    VectorType direction = {
        dx / (dx + dy + dz),
        dy / (dx + dy + dz),
        dz / (dx + dy + dz)
    };

    Scalar epsilon = testEpsilon<Scalar>();
    std::vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise,
                                                     _bUnoriented);
    }

    Ponca::KdTreeDense<DataPoint> kdtree(vectorPoints); // TODO : pass this to the device
    std::cout << "nodes size : " << kdtree.nodes().size() << std::endl;
    // TODO : Try to send the Node to the GPU, or make the kdtree shared memory
    typedef typename Ponca::KdTreeDefaultTraits<DataPoint>::NodeType NodeType;
    NodeType* nodesDevice = nullptr;
    const size_t nodeVectorBufferSize = kdtree.nodes().size() * sizeof(NodeType);
    cudaMalloc(&nodesDevice, nodeVectorBufferSize);
    cudaMemcpy(nodesDevice, kdtree.nodes().data(), nodeVectorBufferSize, cudaMemcpyHostToDevice);

    auto scalarBufferSize = nbPoints*sizeof(Scalar);
    auto vectorBufferSize = scalarBufferSize*Dim;

    // Convert point vector to flattened arrays
    Scalar positions[nbPoints*Dim];
    Scalar normals  [nbPoints*Dim];
    pointsToFlattenedArray<DataPoint>(vectorPoints, positions, normals);

    // Send inputs to the device
    Scalar* positionsDevice;
    Scalar* normalsDevice;
    cudaMalloc(&positionsDevice, vectorBufferSize);
    cudaMalloc(&normalsDevice  , vectorBufferSize);
    cudaMemcpy( positionsDevice, positions, vectorBufferSize, cudaMemcpyHostToDevice);
    cudaMemcpy( normalsDevice  , normals  , vectorBufferSize, cudaMemcpyHostToDevice);

    // Prepare output buffers
    auto *potentialResults = new Scalar[nbPoints];
    auto *gradientResults  = new Scalar[nbPoints*Dim];
    Scalar* potentialResultsDevice;
    Scalar* gradientResultsDevice;
    cudaMalloc(&potentialResultsDevice, scalarBufferSize);
    cudaMalloc(&gradientResultsDevice , vectorBufferSize);

    int blockSize = 128;
    // The grid size needed, based on input size
    int gridSize = (nbPoints + blockSize - 1) / blockSize;

    // Computes the kernel
    fitPotentialKernel<DataPoint, MeanFitSmooth><<<gridSize, blockSize>>>(positionsDevice, normalsDevice, nodesDevice, analysisScale, nbPoints, potentialResultsDevice, gradientResultsDevice);

    cudaDeviceSynchronize(); // Wait for the results

    // Fetch results (Device to Host)
    cudaMemcpy(potentialResults, potentialResultsDevice, scalarBufferSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(gradientResults , gradientResultsDevice , vectorBufferSize, cudaMemcpyDeviceToHost);

    // Free CUDA memory
    cudaFree(positionsDevice);
    cudaFree(normalsDevice);
    cudaFree(potentialResultsDevice);
    cudaFree(gradientResultsDevice);

    for (int j = 0; j < nbPoints; ++j) {
        VectorType primGrad = extractVectorFromFlattenedArray<DataPoint>(j, gradientResults);
        std::cout << "j:" << j << ", potential:"<< potentialResults[j] << ", ";
        // std::cout << "primitiveGradient:" << primGrad.transpose() << " ; " << std::endl;
        VERIFY(std::abs(potentialResults[j]) <= epsilon);
        // VERIFY(Scalar(1.) - std::abs(primGrad.dot(direction)) <= epsilon);
    }

    delete[] potentialResults;
    delete[] gradientResults;
}


__host__ int main(int argc, char** argv) {
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    std::cout << "Test plane fitting on cuda..." << std::endl;
    testPlaneCuda<float, 3>();
}
