#pragma once

#include "baseType.h"
#include <utility>

#include "differentialQuantities.hpp"
#include "Ponca/src/Fitting/enums.h"
#include "Ponca/src/Fitting/build_Estimator2/estimatorFactory.h"
#include "Ponca/src/SpatialPartitioning/KdTree/kdTree.h"

namespace Ponca::Estimators {
    /// Generic processing function: traverse point cloud, compute fitting, and use functor to process fitting output
    /// \note Functor is called only if fit is stable
    template<typename FitT, typename PointContainer, typename Functor>
    void processPointCloud(
        const PointContainer& pointCloud,
        const Functor f, // Functor called on each completion of the fitting process, for each points in the point cloud
        const typename FitT::Scalar scale,
        const int nbMLSIter = 1
    ) {
        #pragma omp parallel for
        for (int i = 0; i < pointCloud.size(); ++i) {
            FitT fit = computeFitForSinglePoint(i, pointCloud, scale, f, nbMLSIter);
        }
    }

    /// Generic processing function: traverse point cloud, compute fitting, and use functor to process fitting output
    /// \note Functor is called only if fit is stable
    template<typename FitT, typename Functor>
    void processTree(
        const KdTree<PPAdapter>& tree,
        const Functor f, // Functor called on each completion of the fitting process, for each points in the point cloud
        const typename FitT::Scalar scale,
        const int nbMLSIter = 1
    ) {
        #pragma omp parallel for
        for (int i = 0; i < tree.samples().size(); ++i) {
            FitT fit = computeFitForSinglePoint(i, tree.points(), scale, f, nbMLSIter, tree.range_neighbors(i, scale));
        }
    }

    /// Generic processing function: compute fitting on a PointContainer : use functor to process fitting output
    /// \note Functor is called only if fit is stable
    template<typename FitT, typename PointContainer, typename Functor, typename PointIterator>
    Ponca::FIT_RESULT computeFitForSinglePoint(
        const int indexQuery,
        const PointContainer& pointCloud,
        const typename FitT::Scalar scale,
        const Functor f, // Function called when the fitting is completed and stable
        const int nbMLSIter = 1,
        const PointIterator& preComputedNeighborhoodIterator = nullptr, // Typically, set this to a kdtree or knnGraph range_neighbors(p, scale)
        FitT fit = FitT()
    ) {
        Ponca::FIT_RESULT res = UNDEFINED;
        typename FitT::VectorType positionQuery = pointCloud[indexQuery].pos();

        for( int mm = 0; mm < nbMLSIter; ++mm) {
            fit.setWeightFunc({ positionQuery, scale });

            if (preComputedNeighborhoodIterator != nullptr)
                res = fit.computeWithIds(preComputedNeighborhoodIterator);
            else
                res = fit.compute(pointCloud);

            if (res != Ponca::FIT_RESULT::STABLE) {
                std::cerr << "Warning: fit " << indexQuery << " is not stable" << std::endl;
                return res;
            }
            positionQuery = fit.project( positionQuery );
        }

        f(indexQuery, fit, positionQuery);
        return res;
    }

    /// Generic processing function: compute fitting on a PointContainer :
    /// \note Functor is called only if fit is stable. Takes in parameter f(i, fit, positionQuery)
    template<typename PointContainer, typename Functor, typename WeightFunc>
    void estimateDifferentialQuantities(
        FitType name,
        const int indexQuery,
        const PointContainer& pointCloud,
        const PPAdapter::Scalar scale,
        const Functor f,
        const int nbMLSIter = 1
    ) {
        EstimatorFactory< PPAdapter, WeightFunc > factory = getEstimatorFactory<PPAdapter, WeightFunc>();
        BasketBase< PPAdapter, WeightFunc > fit = factory->getEstimator(name);

        if (!fit)
            throw std::runtime_error("Invalid fit");

        computeFitForSinglePoint(indexQuery, pointCloud, scale, f, nbMLSIter, fit);
    }

}
