# Defines one Fitting PCH configuration and attaches it to the Fitting library
add_library(PoncaFittingPCH INTERFACE)

target_precompile_headers(PoncaFittingPCH INTERFACE
        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPCH.h
)

#target_include_directories(PoncaFittingPCH PUBLIC
#        "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
#        "$<INSTALL_INTERFACE:include/>"
#)
target_include_directories(PoncaFittingPCH INTERFACE
        ${PONCA_src_ROOT}
)

#target_compile_options(PoncaFittingPCH INTERFACE
#        $<$<CXX_COMPILER_ID:GNU,Clang>:-include ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPCH.h>
#)

target_link_libraries(PoncaFittingPCH INTERFACE Ponca::Fitting)

add_library(Ponca::FittingPCH ALIAS PoncaFittingPCH)
