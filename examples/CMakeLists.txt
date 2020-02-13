OPTION(PONCA_CompilePCLExample "CompilePCLExample" OFF)

add_custom_target(ponca_examples)

add_subdirectory(cpp)
add_subdirectory(cuda)
#add_subdirectory(python)
#
#find_package(OpenGL)
#if(OPENGL_FOUND)
#    add_subdirectory(glviewer)
#endif()
#
#if(CompilePCLExample)
#    add_subdirectory(pcl)
#endif()

