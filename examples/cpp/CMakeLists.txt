project(Ponca_Examples_cpp LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(ponca_basic_SRCS
    ponca_basic_cpu.cpp
)
add_executable(ponca_basic_cpu ${ponca_basic_SRCS})
add_dependencies(ponca-examples ponca_basic_cpu)
target_link_libraries(ponca_basic_cpu PUBLIC Ponca::Precompiled Eigen3::Eigen)

set(ponca_fit_plane_SRCS
    ponca_fit_plane.cpp
)
add_executable(ponca_fit_plane ${ponca_fit_plane_SRCS})
add_dependencies(ponca-examples ponca_fit_plane)
target_link_libraries(ponca_fit_plane PUBLIC Ponca::Precompiled Eigen3::Eigen)

set(ponca_binding_SRCS
    ponca_binding.cpp
)
add_executable(ponca_binding ${ponca_binding_SRCS})
add_dependencies(ponca-examples ponca_binding)
target_link_libraries(ponca_binding PUBLIC Ponca::Precompiled Eigen3::Eigen)

set(ponca_fit_line_SRCS
    ponca_fit_line.cpp
)
add_executable(ponca_fit_line ${ponca_fit_line_SRCS})
add_dependencies(ponca-examples ponca_fit_line)
target_link_libraries(ponca_fit_line PUBLIC Ponca::Precompiled Eigen3::Eigen)

set(ponca_neighbor_search_SRCS
        ponca_neighbor_search.cpp
)
add_executable(ponca_neighbor_search ${ponca_neighbor_search_SRCS})
add_dependencies(ponca-examples ponca_neighbor_search)
target_link_libraries(ponca_neighbor_search PUBLIC Ponca::Precompiled Eigen3::Eigen)

set(ponca_customize_kdtree_SRCS
        ponca_customize_kdtree.cpp
)
add_executable(ponca_customize_kdtree ${ponca_customize_kdtree_SRCS})
add_dependencies(ponca-examples ponca_customize_kdtree)
target_link_libraries(ponca_customize_kdtree PUBLIC Ponca::Precompiled Eigen3::Eigen)

add_subdirectory(pcl)
add_subdirectory(nanoflann)
