project (Ponca-src)

set(src_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Ponca")

set(ponca_core_INCLUDE
    "Ponca/Core/defines.h"
    "Ponca/Core/algebraicSphere.h"
    "Ponca/Core/algebraicSphere.hpp"
    "Ponca/Core/basket.h"
    "Ponca/Core/covariancePlaneFit.h"
    "Ponca/Core/covariancePlaneFit.hpp"
    "Ponca/Core/curvatureEstimation.h"
    "Ponca/Core/curvatureEstimation.hpp"
    "Ponca/Core/curvature.h"
    "Ponca/Core/curvature.hpp"
    "Ponca/Core/enums.h"
    "Ponca/Core/gls.h"
    "Ponca/Core/gls.hpp"
    "Ponca/Core/meanPlaneFit.h"
    "Ponca/Core/meanPlaneFit.hpp"
    "Ponca/Core/mlsSphereFitDer.h"
    "Ponca/Core/mlsSphereFitDer.hpp"
    "Ponca/Core/mongePatch.h"
    "Ponca/Core/mongePatch.hpp"
    "Ponca/Core/orientedSphereFit.h"
    "Ponca/Core/orientedSphereFit.hpp"
    "Ponca/Core/plane.h"
    "Ponca/Core/primitive.h"
    "Ponca/Core/sphereFit.h"
    "Ponca/Core/sphereFit.hpp"
    "Ponca/Core/unorientedSphereFit.h"
    "Ponca/Core/unorientedSphereFit.hpp"
    "Ponca/Core/weightFunc.h"
    "Ponca/Core/weightFunc.hpp"
    "Ponca/Core/weightKernel.h")

add_library(Ponca INTERFACE)

target_include_directories(Ponca INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
)

set_target_properties(Ponca PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_11
)

target_link_libraries(Ponca PUBLIC INTERFACE Eigen3::Eigen)

install(TARGETS Ponca
    EXPORT PoncaTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

add_library(Ponca::Ponca ALIAS Ponca)

install(FILES "ponca.h"
        DESTINATION ${Ponca_INSTALL_INC_DIR}/)
install(FILES "Ponca/core.h"
        DESTINATION ${Ponca_INSTALL_INC_DIR}/Ponca)
install(FILES ${ponca_core_INCLUDE}
        DESTINATION ${Ponca_INSTALL_INC_DIR}/Ponca/Core)

#############################################
# HACK: have the files showing in the IDE, under the name 'ponca-src'
add_custom_target(ponca-src SOURCES "ponca.h" "Ponca/core.h" ${ponca_core_INCLUDE})
