set(CMAKE_INCLUDE_CURRENT_DIR ON)

################################################################################
# Grenaille                                                                    #
################################################################################

set(grenaille_basic_SRCS 
    grenaille_basic_cpu.cpp 
)
add_executable(grenaille_basic_cpu ${grenaille_basic_SRCS})
add_dependencies(grenaille_examples grenaille_basic_cpu)

set(grenaille_fit_plane_SRCS 
    grenaille_fit_plane.cpp 
)
add_executable(grenaille_fit_plane ${grenaille_fit_plane_SRCS})
add_dependencies(grenaille_examples grenaille_fit_plane)

set(grenaille_binding_SRCS
    grenaille_binding.cpp
)
add_executable(grenaille_binding ${grenaille_binding_SRCS})
add_dependencies(grenaille_examples grenaille_binding)

set(grenaille_mesh_curvature_SRCS
    grenaille_mesh_curvature.cpp
)
add_executable(grenaille_mesh_curvature ${grenaille_mesh_curvature_SRCS})
find_package(OpenMP)
if(OPENMP_FOUND)
    target_compile_options(grenaille_mesh_curvature PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(grenaille_mesh_curvature ${OpenMP_CXX_FLAGS})
endif()


add_dependencies(grenaille_examples grenaille_mesh_curvature)

################################################################################
# CUDA                                                                         #
################################################################################

find_package(CUDA)
if(CUDA_FOUND)
    # Dep: FreeImagePlus
    find_package(FREEIMAGEPLUS)
    if(FREEIMAGEPLUS_FOUND)
        include_directories(${FREEIMAGEPLUS_INCLUDE_PATH})

        # Pass options to NVCC
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3)

        # For compilation ...
        # Specify target & source files to compile it from
        cuda_add_executable( ssgls ssgls.cu )

        # For linking ...
        # Specify target & libraries to link it with
        target_link_libraries( ssgls ${FREEIMAGEPLUS_LIBRARIES} )

        # Add images
        add_custom_command(TARGET ssgls POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy_directory
                               ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data )
    endif(FREEIMAGEPLUS_FOUND)
endif(CUDA_FOUND)
