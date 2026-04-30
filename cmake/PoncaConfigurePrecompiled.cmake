set(ponca_Precompiled_PRECOMPILED
    "${PONCA_src_ROOT}/Ponca/src/Precompiled/precompiled.h"
)

set(ponca_Precompiled_INCLUDE
    "${PONCA_src_ROOT}/Ponca"
)

set(ponca_Precompiled_SRC
    "${PONCA_src_ROOT}/Ponca/Instantiate"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate.cpp"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/entrypoint.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/types.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/dims.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/points.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/filters.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/baskets.h"
    "${PONCA_src_ROOT}/Ponca/src/Common/instantiate/instantiate_list/basketsdiff.h"
    )

add_library(Instantiate SHARED)
add_library(Ponca::Instantiate ALIAS Instantiate)
add_dependencies(Instantiate Fitting Common SpatialPartitioning)


target_include_directories(Instantiate PUBLIC
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )

target_sources(Instantiate PUBLIC
    "$<BUILD_INTERFACE:${ponca_Precompiled_SRC}>"
    "$<INSTALL_INTERFACE:>"
    )

target_link_libraries(Instantiate PUBLIC Eigen3::Eigen)
target_compile_definitions(Instantiate PRIVATE _PONCA_COMPILE_DEFINITION)


# Instantiation options
option(PONCA_INSTANTIATE_ALL "Instantiate everything defined in Ponca/Instantiate. This options takes priority over all other options for instantiation." ON)
if (${PONCA_INSTANTIATE_ALL})
    target_compile_definitions(Instantiate PUBLIC _PONCA_INSTANTIATE_ALL)
    # MSVC complains when too much symbols are within a library... 
    target_compile_options(Instantiate PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
    )    
endif()

set(ponca_Instantiate_OPTIONS
    PONCA_INSTANTIATE_DOUBLE
    PONCA_INSTANTIATE_FLOAT
    PONCA_INSTANTIATE_2D
    PONCA_INSTANTIATE_3D
    PONCA_INSTANTIATE_POINTPOSITION
    PONCA_INSTANTIATE_POINTPOSITIONNORMAL
    PONCA_INSTANTIATE_SMOOTHWEIGHT
    PONCA_INSTANTIATE_CONSTANTWEIGHT
    PONCA_INSTANTIATE_NOWEIGHT
    PONCA_INSTANTIATE_SPACEDER
    PONCA_INSTANTIATE_SCALESPACEDER
    )
foreach(param IN LISTS ponca_Instantiate_OPTIONS)
    option(${param} "" OFF)
    if (${param})
        target_compile_definitions(Instantiate PUBLIC _${param})
    endif()
endforeach()



set_target_properties(Instantiate PROPERTIES
  LINKER_LANGUAGE CXX
)


install(TARGETS Instantiate
    EXPORT InstantiateTargets
    LIBRARY DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION  ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCDIR}
)

install(EXPORT InstantiateTargets
  FILE PoncaTargets-Instantiate.cmake
  NAMESPACE Ponca::
  DESTINATION lib/cmake/Ponca
)

export(EXPORT CommonTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/PoncaTargets-Instantiate.cmake
  NAMESPACE Ponca::
)