#add_library(FittingPCH STATIC
#        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPCH.cpp
#)
#
#target_include_directories(FittingPCH PUBLIC
#        "$<BUILD_INTERFACE:${PONCA_src_ROOT}>"
#        "$<INSTALL_INTERFACE:include/>"
#)
#
#target_link_libraries(FittingPCH PUBLIC Fitting)
#
#target_precompile_headers(FittingPCH PUBLIC
#        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingPCH.h
#)
#
#add_library(Ponca::FittingPCH ALIAS FittingPCH)

#
#add_library(Ponca::FittingPCH INTERFACE)
#
#target_precompile_headers(Ponca::FittingPCH INTERFACE
#        ${PONCA_src_ROOT}/Ponca/precompiled/Fitting/fittingETI.h
#)
