/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <pcl/features/feature.h>
#include <Ponca/src/Common/pointTypes.h>

/*!
 * \file examples/cpp/pcl/pcl_wrapper.h
 * \brief Class used for fitting method.
 * Using this approach, ones can use the ponca library with already existing
 * data-structures and without any data-duplication.
 */

using GlsPoint = Ponca::PointPositionNormal<float, 3>;

namespace pcl
{
    /** \brief GlsCurvature estimates local surface curvatures at each 3D point with oriented normal using the Ponca
     * library method.
     *
     * \author Gautier Ciaudo
     */
    template <typename PointInT, typename PointOutT>
    class GlsCurvature : public Feature<PointInT, PointOutT>
    {
        using Feature<PointInT, PointOutT>::feature_name_;
        using Feature<PointInT, PointOutT>::getClassName;
        using Feature<PointInT, PointOutT>::indices_;
        using Feature<PointInT, PointOutT>::input_;
        using Feature<PointInT, PointOutT>::surface_;
        using Feature<PointInT, PointOutT>::k_;
        using Feature<PointInT, PointOutT>::search_radius_;
        using Feature<PointInT, PointOutT>::search_parameter_;

        using PointCloudOut      = typename Feature<PointInT, PointOutT>::PointCloudOut;
        using PointCloudConstPtr = typename Feature<PointInT, PointOutT>::PointCloudConstPtr;

    public:
        /** \brief Empty constructor. */
        GlsCurvature() { feature_name_ = "GlsCurvature"; }

        /** \brief Compute the Growing Least-Squares sphere fit for a given set of points, using their indices,
         * and return the estimated surface curvature.
         * \param cloud the input point cloud
         * \param p_idx The index of the point of evaluation
         * \param indices the point cloud indices that need to be used
         * \param curvature the estimated surface curvature
         */
        void computeCurvature(const PointCloud<PointInT>& cloud, int p_idx, const std::vector<int>& indices,
                              float& curvature);

    protected:
        /** \brief Estimate curvatures for all points given in <setInputCloud (), setIndices ()> using the surface in
         * setSearchSurface () and the spatial locator in setSearchMethod ()
         * \note In situations where not enough neighbors are found, curvature values are set to qNan.
         * \param[out] output the point cloud model dataset that contains surface curvatures
         */
        void computeFeature(PointCloudOut& output) override;

    private:
        /** \brief Make the computeFeature (&Eigen::MatrixXf); inaccessible from outside the class
         * \param[out] output the output point cloud
         */
        static void computeFeatureEigen(PointCloud<Eigen::MatrixXf>& output) {}
    };
} // namespace pcl
