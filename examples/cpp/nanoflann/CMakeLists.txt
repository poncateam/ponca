project(Ponca_Example_Nanoflann LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(ponca_nanoflann ponca_nanoflann.cpp)
target_include_directories(ponca_nanoflann PRIVATE ${PONCA_src_ROOT})
add_dependencies(ponca-examples ponca_nanoflann)
ponca_handle_eigen_dependency(ponca_nanoflann)
