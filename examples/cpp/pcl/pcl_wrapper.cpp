/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <pcl/point_types.h>
#include <pcl/impl/instantiate.hpp>

#include "pcl_wrapper.h"
#include "pcl_wrapper.hpp"

// Instantiations of specific point types
PCL_INSTANTIATE_PRODUCT(GlsCurvature, ((pcl::PointNormal)(pcl::PointXYZRGBNormal)(pcl::PointXYZINormal))((pcl::PointNormal)(pcl::PointXYZRGBNormal)(pcl::PointXYZINormal)));
