project(Ponca_Examples_python LANGUAGES CXX)

include(CheckLanguage)
check_language(CUDA)

################################################################################
# CUDA                                                                         #
################################################################################
if(CMAKE_CUDA_COMPILER)
    set (CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    enable_language(CUDA)

    if(POLICY CMP0104)
        cmake_policy(SET CMP0104 NEW)
    endif()

    if(NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 20)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()
    set(SSGLS_CU ${CMAKE_CURRENT_SOURCE_DIR}/ssgls.cu)

    # Generate PTX
    add_library(ssgls_ptx OBJECT ${SSGLS_CU})
    set_property(TARGET ssgls_ptx PROPERTY CUDA_ARCHITECTURES OFF)

    # Include directories
    target_include_directories(ssgls_ptx PRIVATE ${PONCA_src_ROOT})
    target_link_libraries(ssgls_ptx PUBLIC Eigen3::Eigen)
    set_target_properties(ssgls_ptx PROPERTIES
        CUDA_PTX_COMPILATION ON
    )

    # Custom target to copy outputs
    add_custom_target(ssgls_py ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/ssgls.py ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_OBJECTS:ssgls_ptx> ${CMAKE_CURRENT_BINARY_DIR}/ssgls.ptx
        DEPENDS ssgls_ptx
    )

else()
    message(STATUS "CUDA not found, skipping Ponca_Examples_python")
endif()
