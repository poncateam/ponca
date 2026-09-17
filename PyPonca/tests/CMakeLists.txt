find_package(Python3)

# Test cases generator
add_executable(generate_pytests3d cpp/test3d.cpp)
target_link_libraries(generate_pytests3d PRIVATE Common Fitting)

# Add generators as fixtures 
add_test(
    NAME pyponca_generate_testcases
    COMMAND $<TARGET_FILE:generate_pytests3d> ${CMAKE_BINARY_DIR}/test_cases.json
)
set_tests_properties(pyponca_generate_testcases PROPERTIES FIXTURES_SETUP pyponca_fixture LABELS "pyponca")

# Add python tests
add_test(
    NAME pyponca_test_api
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_api.py
)
set_tests_properties(pyponca_test_api PROPERTIES FIXTURES_SETUP pyponca_fixture LABELS "pyponca")


add_test(
    NAME pyponca_test_accuracy
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_with_env.py 
        --env-name TEST_CASES --env-value ${CMAKE_BINARY_DIR}/test_cases.json ${CMAKE_CURRENT_SOURCE_DIR}/test_accuracy.py
)
set_tests_properties(pyponca_test_accuracy PROPERTIES FIXTURES_SETUP pyponca_fixture LABELS "pyponca")

