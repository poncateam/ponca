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

#############################################
set(gls_tau_SRCS 
    gls_tau.cpp 
)
add_executable(gls_tau ${gls_tau_SRCS})
add_dependencies(buildtests gls_tau)
add_test(gls_tau gls_tau)

#############################################
set(gls_sphere_der_SRCS 
    gls_sphere_der.cpp 
)
add_executable(gls_sphere_der ${gls_sphere_der_SRCS})
add_dependencies(buildtests gls_sphere_der)
add_test(gls_sphere_der gls_sphere_der)

#############################################
set(gls_paraboloid_der_SRCS 
    gls_paraboloid_der.cpp
)
add_executable(gls_paraboloid_der ${gls_paraboloid_der_SRCS})
add_dependencies(buildtests gls_paraboloid_der)
add_test(gls_paraboloid_der gls_paraboloid_der)

#############################################
set(fit_unoriented_SRCS 
    fit_unoriented.cpp
)
add_executable(fit_unoriented ${fit_unoriented_SRCS})
add_dependencies(buildtests fit_unoriented)
add_test(fit_unoriented fit_unoriented)

#############################################
set(gls_compare_SRCS 
    gls_compare.cpp
)
add_executable(gls_compare ${gls_compare_SRCS})
add_dependencies(buildtests gls_compare)
add_test(gls_compare gls_compare)

#############################################
set(plane_primitive_SRCS
    plane_primitive.cpp
)
add_executable(plane_primitive ${plane_primitive_SRCS})
add_dependencies(buildtests plane_primitive)
add_test(plane_primitive plane_primitive)

#############################################
set(fit_plane_SRCS
    fit_plane.cpp
)
add_executable(fit_plane ${fit_plane_SRCS})
add_dependencies(buildtests fit_plane)
add_test(fit_plane fit_plane)

