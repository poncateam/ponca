OPTION(CompilePCLExample "CompilePCLExample" OFF)

add_custom_target(grenaille_examples)

add_subdirectory(cpp)

if(CompilePCLExample)
	add_subdirectory(pcl)
endif()

