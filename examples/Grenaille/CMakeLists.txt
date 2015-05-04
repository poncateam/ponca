OPTION(CompilePCLExample "CompilePCLExample" OFF)

add_subdirectory(cpp)

if(CompilePCLExample)
	add_subdirectory(pcl)
endif()

