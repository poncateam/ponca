set(CMAKE_INCLUDE_CURRENT_DIR ON)

################################################################################
# Grenaille                                                                    #
################################################################################

set(grenaille_basic_SRCS 
    grenaille_basic_cpu.cpp 
)
add_executable(grenaille_basic_cpu ${grenaille_basic_SRCS})
add_dependencies(examples grenaille_basic_cpu)

set(grenaille_fit_plane_SRCS 
    grenaille_fit_plane.cpp 
)
add_executable(grenaille_fit_plane ${grenaille_fit_plane_SRCS})
add_dependencies(examples grenaille_fit_plane)

set(grenaille_binding_SRCS 
    grenaille_binding.cpp 
)
add_executable(grenaille_binding ${grenaille_binding_SRCS})
add_dependencies(examples grenaille_binding)
