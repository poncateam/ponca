set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(PCL 1.3 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

set(pcl_wrapper_SRCS
	main.cpp
    pcl_wrapper.cpp
	pcl_wrapper.hpp
	pcl_wrapper.h
)
add_executable(pcl_wrapper ${pcl_wrapper_SRCS})

target_link_libraries(pcl_wrapper ${PCL_LIBRARIES})
add_dependencies(grenaille_examples pcl_wrapper)

