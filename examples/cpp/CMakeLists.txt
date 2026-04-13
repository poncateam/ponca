project(Ponca_Examples_cpp LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(ponca_basic_SRCS
    ponca_basic_cpu.cpp
)
add_executable(ponca_basic_cpu ${ponca_basic_SRCS})
target_include_directories(ponca_basic_cpu PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_basic_cpu)
target_link_libraries(ponca_basic_cpu PUBLIC Eigen3::Eigen)

set(ponca_fit_plane_SRCS
    ponca_fit_plane.cpp
)
add_executable(ponca_fit_plane ${ponca_fit_plane_SRCS})
target_include_directories(ponca_fit_plane PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_fit_plane)
target_link_libraries(ponca_fit_plane PUBLIC Eigen3::Eigen)

set(ponca_binding_SRCS
    ponca_binding.cpp
)
add_executable(ponca_binding ${ponca_binding_SRCS})
target_include_directories(ponca_binding PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_binding)
target_link_libraries(ponca_binding PUBLIC Eigen3::Eigen)

set(ponca_fit_line_SRCS
    ponca_fit_line.cpp
)
add_executable(ponca_fit_line ${ponca_fit_line_SRCS})
target_include_directories(ponca_fit_line PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_fit_line)
target_link_libraries(ponca_fit_line PUBLIC Eigen3::Eigen)

set(ponca_neighbor_search_SRCS
        ponca_neighbor_search.cpp
)
add_executable(ponca_neighbor_search ${ponca_neighbor_search_SRCS})
target_include_directories(ponca_neighbor_search PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_neighbor_search)
target_link_libraries(ponca_neighbor_search PUBLIC Eigen3::Eigen)

set(ponca_customize_kdtree_SRCS
        ponca_customize_kdtree.cpp
)
add_executable(ponca_customize_kdtree ${ponca_customize_kdtree_SRCS})
target_include_directories(ponca_customize_kdtree PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_customize_kdtree)
target_link_libraries(ponca_customize_kdtree PUBLIC Eigen3::Eigen)

add_subdirectory(pcl)
add_subdirectory(nanoflann)
