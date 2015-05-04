#ifndef PCL_FEATURES_GLS_CURVATURE_IMPL_H_
#define PCL_FEATURES_GLS_CURVATURE_IMPL_H_

#include "pcl_wrapper.h"

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointOutT> void
pcl::GlsCurvature<PointInT, PointOutT>::computeFeature(PointCloudOut &output)
{
    // Allocate enough space to hold the results
    // \note This resize is irrelevant for a radiusSearch ().
    std::vector<int> nn_indices (k_);
    std::vector<float> nn_dists (k_);

    output.is_dense = true;
    // Save a few cycles by not checking every point for NaN/Inf values if the cloud is set to dense
    if (input_->is_dense)
    {
        // Iterating over the entire index vector
        for (size_t idx = 0; idx < indices_->size (); ++idx)
        {
            if (this->searchForNeighbors ((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
            {
                output.points[idx].curvature = std::numeric_limits<float>::quiet_NaN ();
                output.is_dense = false;
                continue;
            }

            computeCurvature (*surface_, (*indices_)[idx], nn_indices, output.points[idx].curvature);
        }
    }
    else
    {
        // Iterating over the entire index vector
        for (size_t idx = 0; idx < indices_->size (); ++idx)
        {
            if (!isFinite ((*input_)[(*indices_)[idx]]) ||
                this->searchForNeighbors ((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
            {
                output.points[idx].curvature = std::numeric_limits<float>::quiet_NaN ();
                output.is_dense = false;
                continue;
            }

            computeCurvature (*surface_, (*indices_)[idx], nn_indices, output.points[idx].curvature);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointOutT> void
pcl::GlsCurvature<PointInT, PointOutT>::computeCurvature(const pcl::PointCloud<PointInT> &cloud, int p_idx, const std::vector<int> &indices, float &curvature)
{
    typedef GlsPoint::Scalar Scalar;
    typedef Grenaille::DistWeightFunc<GlsPoint, Grenaille::SmoothWeightKernel<Scalar> > WeightFunc;
    typedef Grenaille::Basket<GlsPoint, WeightFunc, Grenaille::OrientedSphereFit, Grenaille::GLSParam> Fit;

    Fit fit;
    // Set a weighting function instance using the search radius of the tree as scale
    fit.setWeightFunc(WeightFunc(search_radius_));

    // Set the evaluation position
    fit.init(cloud.points[p_idx].getVector3fMap());

    // Iterate over indices and fit the primitive
    // A GlsPoint instance is generated on the fly to bind the positions to the
    // library representation. No copy is done at this step.
    for(size_t idx = 0; idx < indices.size (); ++idx)
    {
        int id = indices[idx];

        fit.addNeighbor(GlsPoint(cloud.points[id].getVector3fMap(), cloud.points[id].getNormalVector3fMap()));
    }

    // Finalize fitting
    fit.finalize();

    // Test if the fitting ended without errors. Set curvature to qNan otherwise.
    if(fit.isStable())
    {
        curvature = fit.kappa();
    }
    else
    {
        curvature = std::numeric_limits<float>::quiet_NaN ();
    }
}

#define PCL_INSTANTIATE_GlsCurvature(T, OutT) template class PCL_EXPORTS pcl::GlsCurvature<T, OutT>;

#endif // PCL_FEATURES_GLS_CURVATURE_IMPL_H_