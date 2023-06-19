project(Ponca_Examples_cpp LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(ponca_basic_SRCS
    ponca_basic_cpu.cpp
)
add_executable(ponca_basic_cpu ${ponca_basic_SRCS})
target_include_directories(ponca_basic_cpu PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_basic_cpu)
ponca_handle_eigen_dependency(ponca_basic_cpu)

set(ponca_fit_plane_SRCS
    ponca_fit_plane.cpp
)
add_executable(ponca_fit_plane ${ponca_fit_plane_SRCS})
target_include_directories(ponca_fit_plane PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_fit_plane)
ponca_handle_eigen_dependency(ponca_fit_plane)

set(ponca_binding_SRCS
    ponca_binding.cpp
)
add_executable(ponca_binding ${ponca_binding_SRCS})
target_include_directories(ponca_binding PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_binding)
ponca_handle_eigen_dependency(ponca_binding)

set(ponca_fit_line_SRCS
    ponca_fit_line.cpp
)
add_executable(ponca_fit_line ${ponca_fit_line_SRCS})
target_include_directories(ponca_fit_line PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_fit_line)
ponca_handle_eigen_dependency(ponca_fit_line)

set(ponca_neighbor_search_SRCS
    ponca_neighbor_search.cpp
)
add_executable(ponca_neighbor_search ${ponca_neighbor_search_SRCS})
target_include_directories(ponca_neighbor_search PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_neighbor_search)
ponca_handle_eigen_dependency(ponca_neighbor_search)

add_subdirectory(pcl)
add_subdirectory(nanoflann)
