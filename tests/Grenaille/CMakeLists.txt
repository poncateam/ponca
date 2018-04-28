set(CMAKE_INCLUDE_CURRENT_DIR ON)

################################################################################
# Grenaille                                                                    #
################################################################################

add_multi_test(algebraicsphere_primitive.cpp)
add_multi_test(deta_orthogonal_derivatives.cpp)
add_multi_test(fit_radius_curvature_center.cpp)
add_multi_test(gls_tau.cpp)
add_multi_test(gls_sphere_der.cpp)
add_multi_test(gls_paraboloid_der.cpp)
add_multi_test(fit_unoriented.cpp)
add_multi_test(gls_compare.cpp)
add_multi_test(plane_primitive.cpp)
add_multi_test(fit_plane.cpp)
add_multi_test(basket.cpp)
add_multi_test(curvature_plane.cpp)
add_multi_test(curvature_sphere.cpp)
add_multi_test(normal_sphere.cpp)
add_multi_test(normal_plane.cpp)
add_multi_test(dpotential_plane.cpp)
add_multi_test(dpotential_sphere.cpp)
add_multi_test(dnormal_plane.cpp)
add_multi_test(dnormal_sphere.cpp)
add_multi_test(scale_space_consistency.cpp)
