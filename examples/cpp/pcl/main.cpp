/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <boost/thread/thread.hpp>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/features/principal_curvatures.h>
#include <pcl/features/normal_3d.h>
#include<pcl/visualization/pcl_plotter.h>
#include <pcl/common/time.h>

#include "pcl_wrapper.h"


typedef struct
{
    double r,g,b;
}Color;

// Colormap function
Color getColor(float value, float min_value, float max_value)
{
    Color c = {1.0, 1.0, 1.0};
    double dv;
    // Unknown values in our kernel
    if(value == 0.)
    {
        return c;
    }
    // Threshold
    if (value < min_value)
    {
        value = min_value;
    }
    if (value > max_value)
    {
        value = max_value;
    }
    // Interval
    dv = max_value - min_value;
    // Seismic color map like
    if(value < (min_value + 0.5 * dv))
    {
        c.r = 2 * (value - min_value) / dv;
        c.g = 2 * (value - min_value) / dv;
        c.b = 1;
    }
    else
    {
        c.b = 2 - 2 * (value - min_value) / dv;
        c.g = 2 - 2 * (value - min_value) / dv;
        c.r = 1;
    }
    return c;
}

int main()
{
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_without_normals(new pcl::PointCloud<pcl::PointXYZ>);

    // load an object into a point cloud
    if (pcl::io::loadPLYFile<pcl::PointXYZ> ("bun_zipper.ply", *cloud_without_normals) == -1)
    {
        PCL_ERROR ("bun_zipper.ply \n");
        return (-1);
    }

    float radius = 0.005;

    // calculate surface normals with a search radius of 0.01
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(cloud_without_normals);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
    ne.setSearchMethod(tree);
    pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
    ne.setRadiusSearch(radius);
    ne.compute(*normals);

    // concatenate normals with positions into another pointcloud
    pcl::concatenateFields(*cloud_without_normals, *normals, *cloud);

    // compute the gls curvature feature
    pcl::GlsCurvature<pcl::PointNormal, pcl::PointNormal> ne2;
    ne2.setInputCloud(cloud);
    pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal> ());
    ne2.setSearchMethod(tree2);
    pcl::PointCloud<pcl::PointNormal>::Ptr output_cloud (new pcl::PointCloud<pcl::PointNormal>);
    ne2.setRadiusSearch( radius );
    {
      pcl::ScopeTime t1 ("Curvature computation using Ponca");
      ne2.compute(*output_cloud);
    }

    // compare timings with PCL curvature estimation
    pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal> pce;
    pce.setInputCloud(cloud_without_normals);
    pce.setInputNormals(normals);
    pce.setSearchMethod(tree);
    pce.setRadiusSearch(radius);
    pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principalCurvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
    {
      pcl::ScopeTime t1 ("Curvature computation using PCL");
      pce.compute(*principalCurvatures);
    }

    // Use another point cloud for the visualizer (use a colormap to represent the curvature)
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr
        cloud_rgb1 (new pcl::PointCloud<pcl::PointXYZRGB>), //ponca version
        cloud_rgb2 (new pcl::PointCloud<pcl::PointXYZRGB>); //pcl version
    std::vector<double> curvBuffer1 (cloud->size()),  // ponca version
                        curvBuffer2 (cloud->size());  // pcl version
    for(size_t i = 0; i < cloud->size(); ++i)
    {
        // \FIXME There is problem with the magnitude of the derivative computed with Ponca.
        float curvature1 = 1.f / 6000.f * output_cloud->points[i].curvature;
        curvBuffer1[i] = curvature1;
        float curvature2 = 1.f / 0.5f * (principalCurvatures->points[i].pc1 + principalCurvatures->points[i].pc2)/2.f;
        curvBuffer2[i] = curvature2;

        Color c1 = getColor(curvature1, -0.5f, 0.5f);
        Color c2 = getColor(-curvature2, -0.5f, 0.5f);

        pcl::PointXYZRGB point;

        point.x = cloud->points[i].x;
        point.y = cloud->points[i].y;
        point.z = cloud->points[i].z;

        //convert into int
        point.r = c1.r * 255;
        point.g = c1.g * 255;
        point.b = c1.b * 255;
        cloud_rgb1->points.push_back (point);


        point.r = c2.r * 255;
        point.g = c2.g * 255;
        point.b = c2.b * 255;
        cloud_rgb2->points.push_back (point);
    }
//    // visualize curvature plots
//    pcl::visualization::PCLPlotter * plotter = new pcl::visualization::PCLPlotter ();
//    plotter->addHistogramData (curvBuffer2, 100, "PCL");
//    plotter->addHistogramData (curvBuffer1, 100, "Ponca");
//    plotter->setShowLegend (true); //show legends

    // visualize curvatures
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("Ponca - PCL Demo"));
    viewer->setBackgroundColor (0.5, 0.5, 0.5);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb1(cloud_rgb1);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_rgb1, rgb1, "Ponca estimation");
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb2(cloud_rgb2);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_rgb2, rgb2, "PCL estimation");

    // make cloud visible directly
    Eigen::Affine3f tr {Eigen::Affine3f::Identity()};

    tr.translate( Eigen::Vector3f{ 0.1, -0.1, 0.4 } );
    viewer->updatePointCloudPose( "Ponca estimation", tr);

    tr.translate( Eigen::Vector3f{ -0.2, 0.0, 0.0 } );
    viewer->updatePointCloudPose( "PCL estimation", tr);

    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "Ponca estimation");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "PCL estimation");
    viewer->initCameraParameters ();

    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
//        plotter->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}
