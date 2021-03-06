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

set(ponca_mesh_curvature_SRCS
    ponca_mesh_curvature.cpp
)

add_subdirectory(pcl)
