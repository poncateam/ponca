## Defines one Fitting PCH configuration and attaches it to the Fitting library
add_library(PoncaFittingPCH INTERFACE)
target_precompile_headers(PoncaFittingPCH INTERFACE
        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPCH.h
)
target_include_directories(PoncaFittingPCH INTERFACE ${PONCA_src_ROOT})
target_link_libraries(PoncaFittingPCH INTERFACE Ponca::Fitting)
add_library(Ponca::FittingPCH ALIAS PoncaFittingPCH)

# Explicit instantiation for the templated class
add_library(PoncaFittingETI STATIC
        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPrecompiled.cpp
)
target_include_directories(PoncaFittingETI PUBLIC ${PONCA_src_ROOT})
target_link_libraries(PoncaFittingETI PUBLIC
        Ponca::Fitting
)
add_library(Ponca::FittingETI ALIAS PoncaFittingETI)
