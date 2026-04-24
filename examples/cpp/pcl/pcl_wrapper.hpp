/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <Ponca/Precompiled>

#include "pcl_wrapper.h"
#include <pcl/common/point_tests.h> // isFinite

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointOutT>
void pcl::GlsCurvature<PointInT, PointOutT>::computeFeature(PointCloudOut& output)
{
    // Allocate enough space to hold the results
    // \note This resize is irrelevant for a radiusSearch ().
    std::vector<int> nn_indices(k_);
    std::vector<float> nn_dists(k_);

    output.is_dense = true;
    // Save a few cycles by not checking every point for NaN/Inf values if the cloud is set to dense
    if (input_->is_dense)
    {
        // Iterating over the entire index vector
        for (size_t idx = 0; idx < indices_->size(); ++idx)
        {
            if (this->searchForNeighbors((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
            {
                output.points[idx].curvature = std::numeric_limits<float>::quiet_NaN();
                output.is_dense              = false;
                continue;
            }

            computeCurvature(*surface_, (*indices_)[idx], nn_indices, output.points[idx].curvature);
        }
    }
    else
    {
        // Iterating over the entire index vector
        for (size_t idx = 0; idx < indices_->size(); ++idx)
        {
            if (!pcl::isFinite((*input_)[(*indices_)[idx]]) ||
                this->searchForNeighbors((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
            {
                output.points[idx].curvature = std::numeric_limits<float>::quiet_NaN();
                output.is_dense              = false;
                continue;
            }

            computeCurvature(*surface_, (*indices_)[idx], nn_indices, output.points[idx].curvature);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointOutT>
void pcl::GlsCurvature<PointInT, PointOutT>::computeCurvature(const pcl::PointCloud<PointInT>& cloud, int p_idx,
                                                              const std::vector<int>& indices, float& curvature)
{
    using Scalar     = GlsPoint::Scalar;
    using WeightFunc = Ponca::DistWeightFunc<GlsPoint, Ponca::SmoothWeightKernel<Scalar>>;
    using FitBasket  = Ponca::Basket<GlsPoint, WeightFunc, Ponca::CovariancePlaneFit>;
    using Fit =
        Ponca::BasketDiff<FitBasket, Ponca::FitScaleSpaceDer, Ponca::CovariancePlaneDer, Ponca::CurvatureEstimatorDer,
                          Ponca::NormalDerivativeWeingartenEstimator, Ponca::WeingartenCurvatureEstimatorDer>;

    Fit fit;
    // Set a weighting function instance using the search radius of the tree as scale
    fit.setNeighborFilter({cloud.points[p_idx].getVector3fMap(), float(search_radius_)});

    fit.init();

    // Iterate over indices and fit the primitive
    // A GlsPoint instance is generated on the fly to bind the positions to the
    // library representation. No copy is done at this step.
    for (int id : indices)
    {
        fit.addNeighbor(GlsPoint(cloud.points[id].getVector3fMap(), cloud.points[id].getNormalVector3fMap()));
    }

    // Finalize fitting
    fit.finalize();

    // Test if the fitting ended without errors. Set curvature to qNan otherwise.
    if (fit.isStable())
    {
        curvature = fit.kMean();
    }
    else
    {
        curvature = std::numeric_limits<float>::quiet_NaN();
    }
}

#define PCL_INSTANTIATE_GlsCurvature(T, OutT) template class PCL_EXPORTS pcl::GlsCurvature<T, OutT>;
