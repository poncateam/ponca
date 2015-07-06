OPTION(CompilePCLExample "CompilePCLExample" OFF)

add_custom_target(grenaille_examples)

add_subdirectory(cpp)
add_subdirectory(python)

if(CompilePCLExample)
	add_subdirectory(pcl)
endif()

