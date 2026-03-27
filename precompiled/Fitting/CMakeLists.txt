# STATIC lib for explicit instantiation of the templated class (makes them reusable with the vtable)
add_library(PoncaFittingPrecompiled STATIC fittingPrecompiled.cpp)
# OPTIONAL : To speed up header parsing speed when recompiling a target
if(PONCA_GENERATE_PCH)
    target_precompile_headers(PoncaFittingPrecompiled PUBLIC fittingPCH.h)
endif()
target_include_directories(PoncaFittingPrecompiled PUBLIC ${PONCA_src_ROOT})
target_link_libraries(PoncaFittingPrecompiled PUBLIC Ponca::Fitting)
add_library(Ponca::FittingPrecompiled ALIAS PoncaFittingPrecompiled)
