project(Ponca_Examples_sycle LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(AdaptiveCpp CONFIG QUIET)

if(AdaptiveCpp_FOUND)
    message(STATUS "Compiling Sycl examples")
    
    add_executable(ponca_basic_sycl ponca_basic_sycl.cpp)
    target_include_directories(ponca_basic_sycl PRIVATE ${PONCA_src_ROOT})
    target_link_libraries(ponca_basic_sycl PUBLIC Eigen3::Eigen)
    add_dependencies(ponca-examples ponca_basic_sycl)

    # 3. Add SYCL compilation support to the target
    add_sycl_to_target(TARGET ponca_basic_sycl SOURCES ponca_basic_sycl.cpp)

    add_executable(ponca_fit_kdtree_sycl ponca_fit_kdtree.cpp)
    target_include_directories(ponca_fit_kdtree_sycl PRIVATE ${PONCA_src_ROOT})
    target_link_libraries(ponca_fit_kdtree_sycl PUBLIC Eigen3::Eigen)
    add_dependencies(ponca-examples ponca_fit_kdtree_sycl)

    # 3. Add SYCL compilation support to the target
    add_sycl_to_target(TARGET ponca_fit_kdtree_sycl SOURCES ponca_fit_kdtree.cpp)
else()
    message(STATUS "Sycl not found. Skipping examples")
endif(AdaptiveCpp_FOUND)
