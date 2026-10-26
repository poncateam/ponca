option(PYPONCA_ENABLE_CUDA "Compile PyPonca with CUDA support when a CUDA compiler is found" ON)

if (PYPONCA_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)

    if( CMAKE_CUDA_COMPILER )
        # nvcc must use the same host compiler as the C++ code (and may default to
        # another one, such as on Ubuntu)
        set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
        if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES "native")
        endif()
        if (NOT DEFINED CMAKE_CUDA_STANDARD)
            set(CMAKE_CUDA_STANDARD 20)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        endif()
        enable_language(CUDA)
    else()
        message(STATUS "CUDA compiler not found, PyPonca is compiled without CUDA support")
    endif()
endif()


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

if ( PYPONCA_ENABLE_CUDA AND CMAKE_CUDA_COMPILER )
    set_source_files_properties(${PYPONCA_SRC} PROPERTIES LANGUAGE CUDA)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)

    # OpenMP::OpenMP_CXX only sets flags for C++ sources, host code compiled by nvcc needs them too
    if(OpenMP_CXX_FOUND)
        target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>)
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Nanobind hides the internal symbols of its code. This causes a lot
    # of warnings (that are harmless), this is a way to disable them. 
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wno-attributes>)
    target_compile_options(_pyponca PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wno-attributes>)
endif()

install(TARGETS _pyponca DESTINATION pyponca)
install(FILES pyponca/__init__.py DESTINATION pyponca)