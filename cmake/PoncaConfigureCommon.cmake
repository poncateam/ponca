set(ponca_Common_INCLUDE
    "${PONCA_src_ROOT}/Ponca/Common"
    "${PONCA_src_ROOT}/Ponca/Ponca"
    "${PONCA_src_ROOT}/Ponca/src/Common/defines.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/Containers/limitedPriorityQueue.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/Containers/stack.h"
    )

add_library(Common INTERFACE)
target_include_directories(Common INTERFACE
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )
target_sources(Common INTERFACE
    "$<BUILD_INTERFACE:${ponca_Fitting_INCLUDE}>"
    "$<INSTALL_INTERFACE:>"
    )


set_target_properties(Common PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_11
)

#target_link_libraries(Fitting PUBLIC INTERFACE Eigen3::Eigen)

install(TARGETS Common
    EXPORT CommonTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT CommonTargets
  FILE PoncaTargets-Common.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake
)

add_library(Ponca::Common ALIAS Common)

#############################################
# HACK: have the files showing in the IDE, under the name 'ponca-src'
#target_sources(Fitting INTERFACE $<BUILD_INTERFACE:${ponca_INCLUDE}> )
add_custom_target(ponca_Common_IDE SOURCES ${ponca_Common_INCLUDE})
