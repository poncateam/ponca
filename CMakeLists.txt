cmake_minimum_required(VERSION 3.10)

# add this options before PROJECT keyword
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(Ponca LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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


# get changeset id
find_package(Git)
if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${Ponca_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_CHANGESET)
endif()

################################################################################
## Dependencies
################################################################################

## Eigen: automatic configuration:
##  1. if EIGEN3_INCLUDE_DIR is set, use it directly
##  2. if any, and version >= 3.3.x, use system version
##  3. otherwise, download (if required) and use submodule
if( NOT EIGEN3_INCLUDE_DIR )
    find_package(Eigen3 QUIET)
    if( (NOT Eigen3_FOUND) OR EIGEN_VERSION_NUMBER VERSION_LESS 3.3 )

        set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/external/eigen")

        if( NOT EXISTS ${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library )
            execute_process(COMMAND git submodule update --init -- ${EIGEN3_INCLUDE_DIR}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        endif( NOT EXISTS ${EIGEN3_INCLUDE_DIR}/signature_of_eigen3_matrix_library )
    endif()
endif( NOT EIGEN3_INCLUDE_DIR )
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "Eigen3 root path: ${EIGEN3_INCLUDE_DIR}")

################################################################################
# Add core code source, so it is visible in dev environements                  #
################################################################################

add_subdirectory( src )

################################################################################
# Generate install target                                                      #
################################################################################
install(EXPORT PoncaTargets
  FILE
    PoncaTargets.cmake
  NAMESPACE
    Ponca::
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake
)


set(config_install_dir "${CMAKE_INSTALL_PREFIX}/lib/cmake/")
set(include_install_dir "${CMAKE_INSTALL_PREFIX}/include")
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

## Configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
#set(targets_export_name "${PROJECT_NAME}Targets")
#et(namespace "${PROJECT_NAME}::")

## Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file(
    "${version_config}"
                VERSION ${Ponca_VERSION}
                COMPATIBILITY AnyNewerVersion
)

## Configure '<PROJECT-NAME>Config.cmake'
## version to be installed
## Use variables:
##   * targets_export_name
##   * PROJECT_NAME
## Used by Config.cmake.in
set(Ponca_INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(Ponca_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set (Ponca_EXPORT_TARGET_DIR \${CMAKE_CURRENT_LIST_DIR})
configure_package_config_file(
    "cmake/Config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_install_dir}"
)

#install(EXPORT Ponca-srcTargets  DESTINATION ${config_install_dir} NAMESPACE gr::)


# Config
#   * <prefix>/lib/cmake/Foo/FooConfig.cmake
#   * <prefix>/lib/cmake/Foo/FooConfigVersion.cmake
install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${config_install_dir}"
)


################################################################################
# add a target to generate API documentation with Doxygen                      #
################################################################################

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set (PONCA_PROJECT_NUMBER ${GIT_CHANGESET})

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
    # Do NOT use list( APPEND ... ) with flags, as adding flags to a non-empty
    # CMAKE_[...]_FLAGS introduces unwanted semicolons.
    set( CMAKE_CXX_FLAGS "-Wall -pedantic -Wextra  -Wdouble-promotion ${CMAKE_CXX_FLAGS}")
endif (CMAKE_COMPILER_IS_GNUCXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#add_subdirectory(examples EXCLUDE_FROM_ALL)

#include(cmake/ConfigureTesting.cmake)
#add_subdirectory(tests EXCLUDE_FROM_ALL)

