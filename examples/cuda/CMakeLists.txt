project(Ponca_Examples_cuda LANGUAGES CXX)

include(CheckLanguage)
check_language(CUDA)

    if( CMAKE_CUDA_COMPILER )

    enable_language(CUDA)

    if(NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 11)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()



    find_package(PNG QUIET)
    if ( PNG_FOUND )
        add_executable(ponca_ssgls "ponca_ssgls.cu")
        target_include_directories(ponca_ssgls PRIVATE ${PONCA_ROOT_DIR})
        target_link_libraries(ponca_ssgls PRIVATE PNG::PNG)
        target_compile_options(ponca_ssgls PRIVATE --expt-relaxed-constexpr)
        add_dependencies(ponca_examples ponca_ssgls)

        # Copy assets
        add_custom_command(TARGET ponca_ssgls POST_BUILD
                           COMMAND ${CMAKE_COMMAND} -E copy_directory
                               ${CMAKE_CURRENT_SOURCE_DIR}/data
                               $<TARGET_FILE_DIR:ponca_ssgls>/data
                           COMMENT "Copying ssgls data to build tree"
                           VERBATIM
                           )
    else()
        message("LibPNG not found, skipping Ponca_ssgls")
    endif ( PNG_FOUND )

else()
    message(INFO "Cuda not found, skipping Ponca_Examples_cuda")
endif( CMAKE_CUDA_COMPILER )


#################################################################################
## CUDA                                                                         #
#################################################################################

#find_package(CUDA QUIET)
#if(CUDA_FOUND)
#    # Dep: FreeImagePlus
#    find_package(FREEIMAGEPLUS QUIET)
#    if(FREEIMAGEPLUS_FOUND)
#        include_directories(${FREEIMAGEPLUS_INCLUDE_PATH})

#        # Pass options to NVCC
#        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3)

#        # For compilation ...
#        # Specify target & source files to compile it from
#        cuda_add_executable( ssgls ssgls.cu )

#        # For linking ...
#        # Specify target & libraries to link it with
#        target_link_libraries( ssgls ${FREEIMAGEPLUS_LIBRARIES} )
#        target_include_directories(ssgls PRIVATE ${PONCA_ROOT_DIR})

#        # Add images
#        add_custom_command(TARGET ssgls POST_BUILD
#                       COMMAND ${CMAKE_COMMAND} -E copy_directory
#                               ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data )
#    endif(FREEIMAGEPLUS_FOUND)
#endif(CUDA_FOUND)