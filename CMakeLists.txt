project(Patate)

cmake_minimum_required(VERSION 2.8)

# guard against in-source builds (source: Eigen)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# guard against bad build-type strings (source: Eigen)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

################################################################################
# Generate Vitelotte's shaders                                                 #
# TODO: find a more appropriate place for this.                                #
################################################################################

set(vitelotte_SHADERS
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/vert_common.glsl
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/geom_common.glsl
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/frag_common.glsl
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/frag_linear.glsl
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/frag_quadratic.glsl
    ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders/frag_wireframe.glsl
)

add_custom_command(
    OUTPUT ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders.hpp
    COMMAND python3 ${Patate_SOURCE_DIR}/bin/shader2cpp
        -l
        -n Vitelotte::VGMeshRendererShaders
        -g _VITELOTTE_UTILS_VG_MESH_RENDERER_SHADERS_
        -o ${Patate_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders.hpp
        ${vitelotte_SHADERS}
    DEPENDS ${vitelotte_SHADERS}
    WORKING_DIRECTORY ${Patate_SOURCE_DIR}
    VERBATIM
)

add_custom_target(vitelottes_vg_mesh_renderer_shaders
    DEPENDS ${PROJECT_SOURCE_DIR}/Patate/Vitelotte/Utils/vgMeshRendererShaders.hpp)

################################################################################
# add a target to generate API documentation with Doxygen                      #
################################################################################

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc-src/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  add_custom_target(doc
     ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doc-src
     COMMAND ${CMAKE_COMMAND} -E copy_directory
       ${CMAKE_CURRENT_SOURCE_DIR}/doc-src/images/
       ${CMAKE_CURRENT_BINARY_DIR}/doc/html/images
     COMMAND ${CMAKE_COMMAND} -E copy
       ${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG
       ${CMAKE_CURRENT_BINARY_DIR}/doc-src/
     COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
     DEPENDS
       ${CMAKE_CURRENT_SOURCE_DIR}/doc-src/Doxyfile.in
       ${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG
     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
     COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

################################################################################
# tests and examples                                                           #
################################################################################

# enable extra warnings
if (CMAKE_COMPILER_IS_GNUCXX)
    list( APPEND CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra  -Wdouble-promotion ${CMAKE_CXX_FLAGS}")
endif (CMAKE_COMPILER_IS_GNUCXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_subdirectory(examples EXCLUDE_FROM_ALL)

include(cmake/ConfigureTesting.cmake)
add_subdirectory(tests EXCLUDE_FROM_ALL)

