cmake_minimum_required(VERSION 3.10)

# add this options before PROJECT keyword
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(Ponca)

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
# Add core code source, so it is visible in dev environements                  #
################################################################################

FILE ( GLOB PONCA_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Ponca/*.h
)

FILE ( GLOB PONCA_IMPL
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Ponca/*.hpp
)

add_library (Ponca OBJECT ${PONCA_HEADERS} ${PONCA_IMPL})
set_target_properties(Ponca PROPERTIES LINKER_LANGUAGE C++)

################################################################################
# Generate install target                                                      #
################################################################################
if (NOT INCLUDE_INSTALL_DIR)
    set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif (NOT INCLUDE_INSTALL_DIR)
install(DIRECTORY Ponca DESTINATION ${INCLUDE_INSTALL_DIR} FILES_MATCHING REGEX ".*\\.(hpp|h)$")
message("Ponca will be installed in ${INCLUDE_INSTALL_DIR}/Ponca")


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

