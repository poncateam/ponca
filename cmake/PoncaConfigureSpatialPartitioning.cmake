set(ponca_SpatialPartitioning_INCLUDE
    "${PONCA_src_ROOT}/Ponca/SpatialPartitioning"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/defines.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/query.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/indexSquaredDistance.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTree.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTree.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTreeNode.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTreeQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeKNearestIndexQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeKNearestIndexQuery.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeKNearestPointQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeKNearestPointQuery.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeNearestIndexQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeNearestIndexQuery.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeNearestPointQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeNearestPointQuery.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeRangeIndexQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeRangeIndexQuery.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeRangePointQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeRangePointQuery.hpp"
    )

add_library(SpatialPartitioning INTERFACE)
target_include_directories(SpatialPartitioning INTERFACE
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )
target_sources(SpatialPartitioning INTERFACE
    "$<BUILD_INTERFACE:${ponca_SpatialPartitioning_INCLUDE}>"
    "$<INSTALL_INTERFACE:>"
    )
add_dependencies(SpatialPartitioning Common)

set_target_properties(SpatialPartitioning PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_11
)

target_link_libraries(SpatialPartitioning PUBLIC INTERFACE Eigen3::Eigen)

install(TARGETS SpatialPartitioning
    EXPORT SpatialPartitioningTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT SpatialPartitioningTargets
  FILE PoncaTargets-SpatialPartitioning.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake
  COMPONENT Common
)

add_library(Ponca::SpatialPartitioning ALIAS SpatialPartitioning)

#############################################
# HACK: have the files showing in the IDE, under the name 'ponca-src'
#target_sources(SpatialPartitioning INTERFACE $<BUILD_INTERFACE:${ponca_INCLUDE}> )
add_custom_target(ponca_SpatialPartitioning_IDE SOURCES ${ponca_SpatialPartitioning_INCLUDE})
