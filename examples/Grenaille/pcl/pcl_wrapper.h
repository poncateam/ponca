#ifndef PCL_FEATURES_GLS_CURVATURE_H_
#define PCL_FEATURES_GLS_CURVATURE_H_

#include <pcl/features/feature.h>

#include "Patate/grenaille.h"
#include "Eigen/Eigen"

/**
   \brief Class used for fitting method.
  
   Using this approach, ones can use the patate library with already existing
   data-structures and without any data-duplication.
 */
class GlsPoint
{
public:
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1>   VectorType;

    MULTIARCH inline GlsPoint(Eigen::Map< const VectorType > pos, Eigen::Map< const VectorType > normal)
        : pos_   (pos), 
          normal_(normal)
    {}

    MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return pos_; }  
    MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return normal_; }

private:
    Eigen::Map< const VectorType > pos_, normal_;
};

namespace pcl
{
    /** \brief GlsCurvature estimates local surface curvatures at each 3D point with oriented normal using the Patate library
    * method.
    *
    *
    * \author Gautier Ciaudo
    * \ingroup features
    */
    template<typename PointInT, typename PointOutT>
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

        typedef typename Feature<PointInT, PointOutT>::PointCloudOut PointCloudOut;
        typedef typename Feature<PointInT, PointOutT>::PointCloudConstPtr PointCloudConstPtr;

    public:
        /** \brief Empty constructor. */
        GlsCurvature()
        {
            feature_name_ = "GlsCurvature";
        }

        /** \brief Compute the Growing Least-Squares sphere fit for a given set of points, using their indices,
        * and return the estimated surface curvature.
        * \param cloud the input point cloud
        * \param indices the point cloud indices that need to be used
        * \param curvature the estimated surface curvature
        */
        void computeCurvature(const pcl::PointCloud<PointInT> &cloud, int p_idx, const std::vector<int> &indices, float &curvature);

    protected:
        /** \brief Estimate curvatures for all points given in <setInputCloud (), setIndices ()> using the surface in
        * setSearchSurface () and the spatial locator in setSearchMethod ()
        * \note In situations where not enough neighbors are found, curvature values are set to qNan.
        * \param output the resultant point cloud model dataset that contains surface curvatures
        */
        void
        computeFeature (PointCloudOut &output);

    private:
        /** \brief Make the computeFeature (&Eigen::MatrixXf); inaccessible from outside the class
          * \param[out] output the output point cloud
          */
        void
        computeFeatureEigen (pcl::PointCloud<Eigen::MatrixXf> &) {}
    };
}
#endif // PCL_FEATURES_GLS_CURVATURE_H_