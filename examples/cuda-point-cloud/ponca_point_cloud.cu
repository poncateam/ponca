
#include "../../tests/common/testUtils.h"
#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/SpatialPartitioning>
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <vector>

//
// template<typename DataPoint>
// typename DataPoint::Scalar generateData(Ponca::KdTree<DataPoint>& tree)
// {
//     typedef typename DataPoint::Scalar Scalar;
//     typedef typename DataPoint::VectorType VectorType;
//
//     //generate sampled sphere
// #ifdef NDEBUG
//     int nbPoints = Eigen::internal::random<int>(500, 1000);
// #else
//     int nbPoints = Eigen::internal::random<int>(100, 200);
// #endif
//
//     Scalar radius = Eigen::internal::random<Scalar>(1., 10.);
//
//     Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
//     Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
//     VectorType center = VectorType::Random() * centerScale;
//
//     std::vector<DataPoint> vectorPoints(nbPoints);
//
// #ifdef NDEBUG
// #pragma omp parallel for
// #endif
//     for(int i = 0; i < int(vectorPoints.size()); ++i)
//     {
//         vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, false, false, false);
//     }
//
//     tree.clear();
//     tree.build(vectorPoints);
//
//     return analysisScale;
// }

template<typename Fit>
__global__ void kernel()
{
    // Fit fit;
    // fit.setWeightFunc({ });
    // fit.init();
}

int main()
{}
