project (Ponca-src)

set(src_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Ponca")

set(ponca_INCLUDE
    "Ponca/algebraicSphere.h"
    "Ponca/algebraicSphere.hpp"
    "Ponca/basket.h"
    "Ponca/covariancePlaneFit.h"
    "Ponca/covariancePlaneFit.hpp"
    "Ponca/curvatureEstimation.h"
    "Ponca/curvatureEstimation.hpp"
    "Ponca/curvature.h"
    "Ponca/curvature.hpp"
    "Ponca/enums.h"
    "Ponca/gls.h"
    "Ponca/gls.hpp"
    "Ponca/meanPlaneFit.h"
    "Ponca/meanPlaneFit.hpp"
    "Ponca/mlsSphereFitDer.h"
    "Ponca/mlsSphereFitDer.hpp"
    "Ponca/mongePatch.h"
    "Ponca/mongePatch.hpp"
    "Ponca/orientedSphereFit.h"
    "Ponca/orientedSphereFit.hpp"
    "Ponca/plane.h"
    "Ponca/primitive.h"
    "Ponca/sphereFit.h"
    "Ponca/sphereFit.hpp"
    "Ponca/unorientedSphereFit.h"
    "Ponca/unorientedSphereFit.hpp"
    "Ponca/weightFunc.h"
    "Ponca/weightFunc.hpp"
    "Ponca/weightKernel.h")

add_library(Ponca INTERFACE)
#target_sources(Ponca INTERFACE
#    $<INSTALL_INTERFACE:${ponca_INCLUDE}>
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${ponca_INCLUDE}>
#    )

target_include_directories(Ponca INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
)

target_link_libraries(Ponca PUBLIC INTERFACE Eigen3::Eigen)

include(GNUInstallDirs)
install(TARGETS Ponca
    EXPORT PoncaTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

add_library(Ponca::Ponca ALIAS Ponca)

install(FILES ${ponca_INCLUDE}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Ponca)

#############################################
#
# HACK: have the files showing in the IDE, under the name 'ponca_src_IDE'
# Source: https://stackoverflow.com/questions/27039019/how-to-have-cmake-show-headers-that-are-not-part-of-any-binary-target-in-the-ide/29214327#29214327
#
add_custom_target(ponca_src_IDE SOURCES ${ponca_INCLUDE})
