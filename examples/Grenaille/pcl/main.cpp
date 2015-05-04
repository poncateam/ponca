#include <boost/thread/thread.hpp>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/features/normal_3d.h>

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

    // calculate surface normals with a search radius of 0.01
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(cloud_without_normals);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
    ne.setSearchMethod(tree);
    pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
    ne.setRadiusSearch(0.01);
    ne.compute(*normals);

    // concatenate normals with positions into another pointcloud
    pcl::concatenateFields(*cloud_without_normals, *normals, *cloud);

    // compute the gls curvature feature
    pcl::GlsCurvature<pcl::PointNormal, pcl::PointNormal> ne2;
    ne2.setInputCloud(cloud);
    pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal> ());
    ne2.setSearchMethod(tree2);
    pcl::PointCloud<pcl::PointNormal>::Ptr output_cloud (new pcl::PointCloud<pcl::PointNormal>);
    ne2.setRadiusSearch(0.02);
    ne2.compute(*output_cloud);

    // Use another point cloud for the visualizer (use a colormap to represent the curvature)
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_rgb (new pcl::PointCloud<pcl::PointXYZRGB>);
    for(int i = 0; i < cloud->size(); ++i)
    {
        float curvature = output_cloud->points[i].curvature;
        
        Color c = getColor(curvature, -80.f, 50.f);

        pcl::PointXYZRGB point;

        point.x = cloud->points[i].x;
        point.y = cloud->points[i].y;
        point.z = cloud->points[i].z;

        //convert into int
        point.r = c.r * 255;
        point.g = c.g * 255;
        point.b = c.b * 255;

        cloud_rgb->points.push_back (point);
    }

    // visualize curvatures
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0.5, 0.5, 0.5);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud_rgb);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_rgb, rgb, "sample cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
    viewer->initCameraParameters ();

    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}