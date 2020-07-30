project(ponca-tests)

include_directories(${CMAKE_SOURCE_DIR})

# Generate split test header file only if it does not yet exist
# in order to prevent a rebuild everytime cmake is configured.
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h "")
  foreach(i RANGE 1 999)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h
      "#ifdef TEST_PART_${i}\n"
      "#define CALL_SUBTEST_${i}(FUNC) CALL_SUBTEST(FUNC)\n"
      "#else\n"
      "#define CALL_SUBTEST_${i}(FUNC)\n"
      "#endif\n\n"
    )
  endforeach()
endif()

# similar to set_target_properties but append the property instead of overwriting it
macro(ponca_add_target_property target prop value)

  get_target_property(previous ${target} ${prop})
  # if the property wasn't previously set, ${previous} is now "previous-NOTFOUND" which cmake allows catching with plain if()
  if(NOT previous)
    set(previous "")
  endif(NOT previous)
  set_target_properties(${target} PROPERTIES ${prop} "${previous} ${value}")
endmacro(ponca_add_target_property)

macro(ponca_configure_compiler_flags target)
  if( CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(PONCA_COVERAGE_TESTING)
      # Add required flags (GCC & LLVM/Clang)
      target_compile_options(${target} PUBLIC
        -O0        # no optimization
        -g         # generate debug info
        --coverage # sets all required flags
        -fprofile-arcs -ftest-coverage
        -DPONCA_COVERAGE_ENABLED
      )
      if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
        target_link_options(${target} PUBLIC --coverage)
      else()
        target_link_libraries(${target} PUBLIC --coverage)
      endif()
    else() # no coverage: enable march
        target_compile_options(${target} PUBLIC -march=native)
    endif()
  endif()
endmacro()

# Macro to add a test
#
# the unique mandatory parameter testname must correspond to a file <testname.cpp>
#
# Depending on the contents of that file, this macro can have 2 behaviors,
# see below.
#
# A. Default behavior
#
# this macro adds an executable <testname> as well as a ctest test named <testname> too.
#
# On platforms with bash simply run:
#   "ctest -V" or "ctest -V -R <testname>"
# On other platform use ctest as usual.
#
# B. Multi-part behavior
#
# If the source file matches the regexp
#    CALL_SUBTEST_[0-9]+|TEST_PART_[0-9]+
# then it is interpreted as a multi-part test.
# The test is split into multiple executables
#    test_<testname>_<N>
# where N runs from 1 to the greatest occurence found in the source file. Each of these
# executables is built passing -DTEST_PART_N. This allows to split large tests
# into smaller executables.
#
# Moreover, targets <testname> are still generated, they
# have the effect of building all the parts of the test.
#
# Again, ctest -R allows to run all matching tests.
macro(add_multi_test filename)

  string(REGEX MATCHALL "(.*)\\.c.*" dummy "${filename}")
  set(testname ${CMAKE_MATCH_1})

  file(READ "${filename}" test_source)
  set(parts 0)
  string(REGEX MATCHALL "CALL_SUBTEST_[0-9]+|TEST_PART_[0-9]+|PONCA_SUFFIXES(;[0-9]+)+" occurences "${test_source}")
  string(REGEX REPLACE "CALL_SUBTEST_|TEST_PART_|PONCA_SUFFIXES" "" suffixes "${occurences}")
  list(REMOVE_DUPLICATES suffixes)
  if(suffixes)
    add_custom_target(${testname})
    add_dependencies(buildtests ${testname})
    foreach(suffix ${suffixes})
      add_executable(${testname}_${suffix} ${filename})
      target_include_directories(${testname}_${suffix} PRIVATE ${PONCA_ROOT_DIR})
      ponca_add_target_property(${testname}_${suffix} COMPILE_FLAGS "-DTEST_PART_${suffix}=1")

      if(OpenMP_CXX_FOUND)
        target_link_libraries(${testname}_${suffix} PUBLIC OpenMP::OpenMP_CXX)
      endif(OpenMP_CXX_FOUND)

      ponca_configure_compiler_flags(${testname}_${suffix})

      add_dependencies(${testname} ${testname}_${suffix})
      add_test(${testname}_${suffix} ${testname}_${suffix})

    endforeach(suffix)
  else(suffixes)
    add_executable(${testname} ${filename})
    target_include_directories(${testname} PRIVATE ${PONCA_ROOT_DIR})

    if(OpenMP_CXX_FOUND)
      target_link_libraries(${testname} PUBLIC OpenMP::OpenMP_CXX)
    endif(OpenMP_CXX_FOUND)

    ponca_configure_compiler_flags(${testname})

    add_dependencies(buildtests ${testname})
    add_test(${testname} ${testname})
  endif(suffixes)
endmacro(add_multi_test)

add_subdirectory(src)


