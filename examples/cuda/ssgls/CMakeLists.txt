find_package(PNG QUIET)
if ( PNG_FOUND )

    add_executable(ponca_ssgls "ponca_ssgls.cu")
    add_dependencies(ponca-examples ponca_ssgls)
    set_property(TARGET ponca_ssgls PROPERTY CUDA_ARCHITECTURES OFF)
    target_compile_options(ponca_ssgls PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    target_link_libraries(ponca_ssgls PRIVATE PNG::PNG)
    target_link_libraries(ponca_ssgls PUBLIC Ponca::Precompiled Eigen3::Eigen)

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
