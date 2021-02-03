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
add_multi_test(fit_monge_patch.cpp)
add_multi_test(basket.cpp)
add_multi_test(projection.cpp)
add_multi_test(kdtree.cpp)
add_multi_test(kdtree_index.cpp)
