set(CMAKE_INCLUDE_CURRENT_DIR ON)


################################################################################
# CUDA                                                                         #
################################################################################

find_package(CUDA)
if(CUDA_FOUND)

    # Pass options to NVCC
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3)

    # For compilation ...
    # Specify target & source files to compile it from
    # cuda_add_executable( ssgls_py ssgls.cu )
    cuda_compile_ptx (ptxfile ssgls.cu)


    add_custom_target(ssgls_py ALL
        DEPENDS ${ptxfile} ssgls.cu
        SOURCES ssgls.cu
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/ssgls.py ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
                ${ptxfile} ${CMAKE_CURRENT_BINARY_DIR}/ssgls.ptx )

endif(CUDA_FOUND)
