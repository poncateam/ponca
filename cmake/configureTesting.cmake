add_custom_target(buildtests)
add_custom_target(check COMMAND "ctest")
add_dependencies(check buildtests)

# check whether /bin/bash exists
find_file(PONCA_BIN_BASH_EXISTS "/bin/bash" PATHS "/" NO_DEFAULT_PATH)

# CMake/Ctest does not allow us to change the build command,
# so we have to workaround by directly editing the generated DartConfiguration.tcl file
# save CMAKE_MAKE_PROGRAM
set(CMAKE_MAKE_PROGRAM_SAVE ${CMAKE_MAKE_PROGRAM})
# and set a fake one
set(CMAKE_MAKE_PROGRAM "@PONCA_MAKECOMMAND_PLACEHOLDER@")

# This call activates testing and generates the DartConfiguration.tcl
# This adds another build target, which is test for Makefile generators,
# or RUN_TESTS for integrated development environments (like Visual Studio)
include(CTest)

set(PONCA_TEST_BUILD_FLAGS "" CACHE STRING "Options passed to the build command of unit tests")

# overwrite default DartConfiguration.tcl
# The worarounds are different for each version of the MSVC IDE
set(PONCA_TEST_TARGET buildtests)
if(MSVC_IDE)
  if(CMAKE_MAKE_PROGRAM_SAVE MATCHES "devenv") # devenv
    set(PONCA_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM_SAVE} Ponca.sln /build Release /project ${PONCA_TEST_TARGET}")
  else() # msbuild
    set(PONCA_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM_SAVE} ${PONCA_TEST_TARGET}.vcxproj /p:Configuration=\${CTEST_CONFIGURATION_TYPE}")
  endif()

  # append the build flags if provided
  if(NOT "${PONCA_TEST_BUILD_FLAGS}" MATCHES "^[ \t]*$")
    set(PONCA_BUILD_COMMAND "${PONCA_BUILD_COMMAND} ${PONCA_TEST_BUILD_FLAGS}")
  endif()

  # apply the dartconfig hack ...
  set(PONCA_MAKECOMMAND_PLACEHOLDER "${PONCA_BUILD_COMMAND}\n#")
else()
  # for make and nmake
  set(PONCA_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM_SAVE} ${PONCA_TEST_TARGET} ${PONCA_TEST_BUILD_FLAGS}")
  set(PONCA_MAKECOMMAND_PLACEHOLDER "${PONCA_BUILD_COMMAND}")
endif()

configure_file(${CMAKE_BINARY_DIR}/DartConfiguration.tcl ${CMAKE_BINARY_DIR}/DartConfiguration.tcl)

# restore default CMAKE_MAKE_PROGRAM
set(CMAKE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM_SAVE})

# un-set temporary variables so that it is like they never existed
unset(CMAKE_MAKE_PROGRAM_SAVE)
unset(PONCA_MAKECOMMAND_PLACEHOLDER)

# enable multi-threading
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("Compiling tests with OpenMP")
endif()

option(PONCA_COVERAGE_TESTING "Enable coverage reporting" OFF)


