cmake_minimum_required(VERSION 3.11)

# add this options before PROJECT keyword
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(Ponca LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(PONCA_ROOT_DIR ${PROJECT_SOURCE_DIR}/src)

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
# Generate install target                                                      #
################################################################################
set(Ponca_INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(Ponca_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(Ponca_EXPORT_TARGET_DIR \${CMAKE_CURRENT_LIST_DIR})

install(EXPORT PoncaTargets
  FILE
    PoncaTargets.cmake
  NAMESPACE
    Ponca::
  DESTINATION
    ${Ponca_INSTALL_LIB_DIR}/cmake
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
configure_package_config_file(
    "cmake/Config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_install_dir}"
)

install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${config_install_dir}"
)


################################################################################
# Add core code source, so it is visible in dev environements                  #
################################################################################

add_subdirectory( src )


################################################################################
# API documentation with Doxygen                                               #
################################################################################
add_subdirectory(doc)

################################################################################
# tests and examples                                                           #
################################################################################
add_subdirectory(examples EXCLUDE_FROM_ALL)

#include(cmake/ConfigureTesting.cmake)
#add_subdirectory(tests EXCLUDE_FROM_ALL)

################################################################################
# Add custom target for Continuous Integration scripts                         #
################################################################################
add_custom_target(ponca-ci SOURCES
                    ${CMAKE_CURRENT_SOURCE_DIR}/.travis.yml)
