project(Ponca_Examples_cuda LANGUAGES CXX)

include(CheckLanguage)
check_language(CUDA)

    if( CMAKE_CUDA_COMPILER )

    set (CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    enable_language(CUDA)

    if(POLICY CMP0104)
        cmake_policy(SET CMP0104 NEW)
    endif()

    if(NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()

    add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int)
    add_definitions(--expt-relaxed-constexpr)

    add_executable(cuda_fit_plane "ponca_fit_plane.cu")
    target_include_directories(cuda_fit_plane PRIVATE ${PONCA_src_ROOT})
    add_dependencies(ponca-examples cuda_fit_plane)
    set_property(TARGET cuda_fit_plane PROPERTY CUDA_ARCHITECTURES OFF)
    ponca_handle_eigen_dependency(cuda_fit_plane)

    add_subdirectory(ssgls)
else()
    message(INFO "Cuda not found, skipping Ponca_Examples_cuda")
endif( CMAKE_CUDA_COMPILER )
