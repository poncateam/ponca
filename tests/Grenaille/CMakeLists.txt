set(CMAKE_INCLUDE_CURRENT_DIR ON)

################################################################################
# Grenaille                                                                    #
################################################################################

#############################################
set(deta_orthogonal_derivatives_SRCS 
    deta_orthogonal_derivatives.cpp 
)
add_executable(deta_orthogonal_derivatives ${deta_orthogonal_derivatives_SRCS})
add_dependencies(buildtests deta_orthogonal_derivatives)
add_test(deta_orthogonal_derivatives deta_orthogonal_derivatives)

#############################################
set(fit_radius_curvature_center_SRCS 
    fit_radius_curvature_center.cpp 
)
add_executable(fit_radius_curvature_center ${fit_radius_curvature_center_SRCS})
add_dependencies(buildtests fit_radius_curvature_center)
add_test(fit_radius_curvature_center fit_radius_curvature_center)
