set(ponca_SpatialPartitioning_INCLUDE
    "${PONCA_src_ROOT}/Ponca/SpatialPartitioning"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/defines.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/query.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/indexSquaredDistance.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTree.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTree.hpp"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/kdTreeTraits.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeQuery.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeKNearestQueries.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeNearestQueries.h"
    "${PONCA_src_ROOT}/Ponca/src/SpatialPartitioning/KdTree/Query/kdTreeRangeQueries.h"
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

if(Eigen3_FOUND)
    message("Compiling with installed Eigen package, enable transitive linking (Version ${Eigen3_VERSION}, path: ${Eigen3_DIR})")
    target_link_libraries(SpatialPartitioning PUBLIC INTERFACE Eigen3::Eigen)
endif()

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
if( ${PONCA_GENERATE_IDE_TARGETS} )
    add_custom_target(ponca_SpatialPartitioning_IDE SOURCES ${ponca_SpatialPartitioning_INCLUDE})
endif()
