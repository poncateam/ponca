#include <pcl/point_types.h>
#include <pcl/impl/instantiate.hpp>

#include "pcl_wrapper.h"
#include "pcl_wrapper.hpp"

// Instantiations of specific point types
PCL_INSTANTIATE_PRODUCT(GlsCurvature, ((pcl::PointNormal)(pcl::PointXYZRGBNormal)(pcl::PointXYZINormal))((pcl::PointNormal)(pcl::PointXYZRGBNormal)(pcl::PointXYZINormal)));