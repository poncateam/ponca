# STATIC lib for explicit instantiation of the templated class (makes them reusable with the vtable)
add_library(PoncaSpatialPartitioningPrecompiled STATIC
        ${PONCA_src_ROOT}/Ponca/precompiled/Common/spatialPartitioningPrecompiled.cpp
)
# OPTIONAL : To speed up header parsing speed when recompiling
if(PONCA_GENERATE_PCH)
    target_precompile_headers(PoncaSpatialPartitioningPrecompiled PUBLIC
            ${PONCA_src_ROOT}/Ponca/precompiled/Common/spatialPartitioningPCH.h
    )
endif()
target_include_directories(PoncaSpatialPartitioningPrecompiled PUBLIC ${PONCA_src_ROOT})
target_link_libraries(PoncaSpatialPartitioningPrecompiled PUBLIC Ponca::SpatialPartitioning)
add_library(Ponca::SpatialPartitioningPrecompiled ALIAS PoncaSpatialPartitioningPrecompiled)
