set(ponca_Fitting_INCLUDE
    "${PONCA_src_ROOT}/Ponca/Fitting"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/defines.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/algebraicSphere.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/algebraicSphere.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/basket.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/covarianceLineFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/covarianceLineFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/covariancePlaneFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/covariancePlaneFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/curvatureEstimation.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/curvatureEstimation.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/curvature.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/curvature.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/enums.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/gls.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/gls.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/meanPlaneFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/meanPlaneFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/mlsSphereFitDer.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/mlsSphereFitDer.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/mongePatch.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/mongePatch.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/orientedSphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/orientedSphereFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/plane.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/primitive.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/sphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/sphereFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/unorientedSphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/unorientedSphereFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/weightFunc.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/weightFunc.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/weightKernel.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/linePrimitive.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/surfacePrimitive.h"
    )

add_library(Fitting INTERFACE)
target_include_directories(Fitting INTERFACE
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )
target_sources(Fitting INTERFACE
    "$<BUILD_INTERFACE:${ponca_Fitting_INCLUDE}>"
    "$<INSTALL_INTERFACE:>"
    )
add_dependencies(Fitting Common)

set_target_properties(Fitting PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_11
)

if(Eigen3_FOUND)
  message("Compiling with installed Eigen package, enable transitive linking")
  target_link_libraries(Fitting PUBLIC INTERFACE Eigen3::Eigen)
endif()

install(TARGETS Fitting
    EXPORT FittingTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT FittingTargets
  FILE PoncaTargets-Fitting.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake
  COMPONENT Common
)

add_library(Ponca::Fitting ALIAS Fitting)

#############################################
# HACK: have the files showing in the IDE, under the name 'ponca-src'
#target_sources(Fitting INTERFACE $<BUILD_INTERFACE:${ponca_INCLUDE}> )
add_custom_target(ponca_Fitting_IDE SOURCES ${ponca_Fitting_INCLUDE})
