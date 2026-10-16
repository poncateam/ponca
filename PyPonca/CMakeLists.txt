# For now cuda compilation is deactivated due to 
# behaviour difference in eigen's solver (presumably)
# include(CheckLanguage)
# check_language(CUDA)
# 
# if( CMAKE_CUDA_COMPILER )
#     enable_language(CUDA)
#     set(CMAKE_CUDA_ARCHITECTURES "native")
# endif()


set(PYPONCA_SRC
    binding/Common/RegisterCommon.cpp
    binding/Fitting/RegisterFitting.cpp
    binding/SpatialPartitioning/RegisterSpatialPartitioning.cpp
    binding/pyponca.cpp
)

nanobind_add_module(_pyponca ${PYPONCA_SRC})
target_link_libraries(_pyponca PUBLIC Common Fitting)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_pyponca PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(_pyponca PUBLIC PYPONCA_WITH_OMP)
endif()

if ( CMAKE_CUDA_COMPILER )
    set_source_files_properties(${PYPONCA_SRC} PROPERTIES LANGUAGE CUDA)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Nanobind hides the internal symbols of its code. This causes a lot
    # of warnings (that are harmless), this is a way to disable them. 
    # TODO: Set the class visibility instead ! 
    target_compile_options(_pyponca PRIVATE -Wno-attributes)
endif()

install(TARGETS _pyponca DESTINATION pyponca)
install(FILES pyponca/__init__.py DESTINATION pyponca)