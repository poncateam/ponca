# STATIC lib for explicit instantiation of the templated class (makes them reusable with the vtable)
add_library(PoncaCommonPrecompiled STATIC
        ${PONCA_src_ROOT}/Ponca/precompiled/Common/commonPrecompiled.cpp
)
# OPTIONAL : To speed up header parsing speed when recompiling
if(PONCA_GENERATE_PCH)
    target_precompile_headers(PoncaCommonPrecompiled PUBLIC
            ${PONCA_src_ROOT}/Ponca/precompiled/Common/commonPCH.h
    )
endif()
target_include_directories(PoncaCommonPrecompiled PUBLIC ${PONCA_src_ROOT})
target_link_libraries(PoncaCommonPrecompiled PUBLIC Ponca::SpatialPartitioning)
add_library(Ponca::CommonPrecompiled ALIAS PoncaCommonPrecompiled)
