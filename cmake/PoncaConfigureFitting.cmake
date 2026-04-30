set(ponca_Fitting_INCLUDE
    "${PONCA_src_ROOT}/Ponca/Fitting"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/defines.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/compute.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/concepts.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/enums.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/CNC/cnc.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/CNC/cnc.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/CNC/cncFormulaEigen.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/Filters/frame.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Filters/weightFilter.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Filters/weightFilter.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Filters/weightKernel.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/ComputeScheme/evaluationScheme.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/ComputeScheme/mlsEvaluationScheme.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/ComputeScheme/project.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/basket.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/gls.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/gls.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/heightField.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/heightField.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/mlsSphereFitDer.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/mlsSphereFitDer.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/weingarten.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/extensions/weingarten.hpp"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/tools/covariance.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/tools/covariance.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/tools/dryFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/tools/mean.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/tools/mean.hpp"
    
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/primitive/algebraicSphere.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/primitive/algebraicSphere.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/primitive/basketUnit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/primitive/plane.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/primitive/linePrimitive.h"

    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/covarianceLineFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/covariancePlaneFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/covariancePlaneFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/meanPlaneFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/mongePatch.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/mongePatch.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/orientedSphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/orientedSphereFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/sphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/sphereFit.hpp"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/unorientedSphereFit.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/Basket/fit/unorientedSphereFit.hpp"
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

target_link_libraries(Fitting PUBLIC INTERFACE Eigen3::Eigen)

install(TARGETS Fitting
    EXPORT FittingTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT FittingTargets
  FILE PoncaTargets-Fitting.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake/Ponca
  COMPONENT Common
)

add_library(Ponca::Fitting ALIAS Fitting)
export(EXPORT FittingTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/PoncaTargets-Fitting.cmake
  NAMESPACE Ponca::
)

#############################################
# HACK: have the files showing in the IDE, under the name 'ponca-src'
#target_sources(Fitting INTERFACE $<BUILD_INTERFACE:${ponca_INCLUDE}> )
if( ${PONCA_GENERATE_IDE_TARGETS} )
  add_custom_target(ponca_Fitting_IDE SOURCES ${ponca_Fitting_INCLUDE})
endif()
