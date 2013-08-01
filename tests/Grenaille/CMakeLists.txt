set(CMAKE_INCLUDE_CURRENT_DIR ON)

################################################################################
# Grenaille                                                                    #
################################################################################

set(deta_orthogonal_derivatives_SRCS 
    deta_orthogonal_derivatives.cpp 
)
add_executable(deta_orthogonal_derivatives ${deta_orthogonal_derivatives_SRCS})
add_dependencies(tests deta_orthogonal_derivatives)
