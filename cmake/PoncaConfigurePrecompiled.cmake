set(ponca_Precompiled_PRECOMPILED
    "${PONCA_src_ROOT}/Ponca/Precompiled"
)

set(ponca_Precompiled_INCLUDE
    "${PONCA_src_ROOT}/Ponca/Ponca"
    "${PONCA_src_ROOT}/Ponca/Ponca"
)

set(ponca_Precompiled_SRC
    "${PONCA_src_ROOT}/Ponca/precompiled/precompiled.cpp"
    )

add_library(Precompiled STATIC)
add_library(Ponca::Precompiled ALIAS Precompiled)

target_include_directories(Precompiled PUBLIC 
    "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
    "$<INSTALL_INTERFACE:include/>"
    )

target_sources(Precompiled PUBLIC
    "$<BUILD_INTERFACE:${ponca_Precompiled_SRC}>"
    "$<INSTALL_INTERFACE:/>"
    )
target_precompile_headers(Precompiled PUBLIC
    "$<BUILD_INTERFACE:${ponca_Precompiled_PRECOMPILED}>"
    )

target_compile_definitions(Precompiled PRIVATE _PONCA_COMPILE_DEFINITION)

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