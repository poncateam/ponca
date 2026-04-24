set(ponca_Precompiled_PRECOMPILED
    "${PONCA_src_ROOT}/Ponca/src/Precompiled/precompiled.h"
)

set(ponca_Precompiled_INCLUDE
    "${PONCA_src_ROOT}/Ponca"
)

set(ponca_Precompiled_SRC
    "${PONCA_src_ROOT}/Ponca/Precompiled"
    "${PONCA_src_ROOT}/Ponca/src/Precompiled/instantiate.h"
    "${PONCA_src_ROOT}/Ponca/src/Precompiled/precompiled.h"
    "${PONCA_src_ROOT}/Ponca/src/Precompiled/precompiled.cpp"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/types.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/dims.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/points.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/instantiate/filters.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/instantiate/baskets.h"
    "${PONCA_src_ROOT}/Ponca/src/Fitting/instantiate/basketsdiff.h"
    )

add_library(Precompiled SHARED)
add_library(Ponca::Precompiled ALIAS Precompiled)
add_dependencies(Precompiled Fitting Common SpatialPartitioning)


target_include_directories(Precompiled PUBLIC
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )

target_sources(Precompiled PUBLIC
    "$<BUILD_INTERFACE:${ponca_Precompiled_SRC}>"
    "$<INSTALL_INTERFACE:>"
    )
target_precompile_headers(Precompiled PUBLIC
    "$<BUILD_INTERFACE:${ponca_Precompiled_PRECOMPILED}>"
    )

target_link_libraries(Precompiled PUBLIC Eigen3::Eigen)
target_compile_definitions(Precompiled PRIVATE _PONCA_COMPILE_DEFINITION)


# Instantiation options
option(PONCA_INSTANTIATE_ALL "Instantiate everything defined in Ponca/Precompiled. This options takes priority over all other options for instantiation." ON)
if (${PONCA_INSTANTIATE_ALL})
    target_compile_definitions(Precompiled PUBLIC _PONCA_INSTANTIATE_ALL)
    # MSVC complains when too much symbols are within a library... 
    target_compile_options(Precompiled PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
    )    
endif()

set(ponca_Precompiled_OPTIONS
    PONCA_INSTANTIATE_DOUBLE
    PONCA_INSTANTIATE_FLOAT
    PONCA_INSTANTIATE_2D
    PONCA_INSTANTIATE_3D
    PONCA_INSTANTIATE_POINTPOSITION
    PONCA_INSTANTIATE_POINTPOSITIONNORMAL
    PONCA_INSTANTIATE_SMOOTHWEIGHT
    PONCA_INSTANTIATE_SPACEDER
    )
foreach(param IN LISTS ponca_Precompiled_OPTIONS)
    option(${param} "" OFF)
    if (${param})
        target_compile_definitions(Precompiled PUBLIC _${param})
    endif()
endforeach()



set_target_properties(Precompiled PROPERTIES
  LINKER_LANGUAGE CXX
)


install(TARGETS Precompiled
    EXPORT PrecompiledTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT PrecompiledTargets
  FILE PoncaTargets-Precompiled.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake/Ponca
)

export(EXPORT CommonTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/PoncaTargets-Precompiled.cmake
  NAMESPACE Ponca::
)