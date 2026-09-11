include(CheckLanguage)
check_language(CUDA)

if( CMAKE_CUDA_COMPILER )
    enable_language(CUDA)
endif()


set(PYPONCA_SRC
    binding/Common/RegisterCommon.cu
    binding/Fitting/RegisterFitting.cu
    binding/SpatialPartitioning/RegisterSpatialPartitioning.cu
    binding/pyponca.cu
)

nanobind_add_module(_pyponca ${PYPONCA_SRC})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_pyponca PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(_pyponca PUBLIC PYPONCA_WITH_OMP)
endif()


if ( CMAKE_CUDA_COMPILER )
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
endif()

# We link here manually with stdc++. Nanobind seems to think this 
# is C-code when linking with our interface module...???
target_link_libraries(_pyponca PUBLIC Common Fitting stdc++)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Nanobind hides the internal symbols of its code. This causes a lot
    # of warnings (that are harmless), this is a way to disable them. 
    # TODO: Set the class visibility instead ! 
    target_compile_options(_pyponca PRIVATE -Wno-attributes)
endif()

install(TARGETS _pyponca DESTINATION pyponca)
install(FILES pyponca/__init__.py DESTINATION pyponca)