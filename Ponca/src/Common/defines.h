/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

/**
 *
 * \defgroup common Common module
 * \brief This modules includes general purpose classes and methods.
 *
 */

////////////////////////////////////////////////////////////////////////////////
// Compatibility types, macros, functions
//
#ifdef __CUDACC__
#    include <cuda.h>
#    define PONCA_MULTIARCH __host__ __device__
#    define PONCA_MULTIARCH_HOST __host__
#else
#    define PONCA_MULTIARCH
#    define PONCA_MULTIARCH_HOST

// GCC: compile with -std=c++0x
#    if defined(__GNUC__) && ((__GNUC__ == 4 && __GNUC_MINOR__ >= 6) || (__GNUC__ >= 5))
#        if defined(nullptr_t) || (__cplusplus > 199711L) || defined(HACK_GCC_ITS_CPP0X)
#            define __CPP0X__
#        endif
#    endif

#endif // ifdef __CUDACC__

#ifdef __CUDACC__
#    define PONCA_MULTIARCH_INCLUDE_STD(FILENAME) "defines.h"
// __device__ version of math function are implicitly defined by cuda and are not inside any
// namespaces. However, other classes (such as numeric_limits) are not. We distinguish both
// cases with the two following macros.
// Note: When the new min supported version of cuda will be 13.0, we will be
// able to merge both maccros due to new addition to libcu++ library.
#    define PONCA_MULTIARCH_CU_STD_FUNC(FUNC) using cuda::std::FUNC;
#    define PONCA_MULTIARCH_STD_MATH(FUNC)
#    define PONCA_MULTIARCH_STD_MATH_NAMESPACE(FUNC) FUNC

// see https://nvidia.github.io/libcudacxx/standard_api.html
#    define PONCA_MULTIARCH_CU_STD_NAMESPACE(FUNC) cuda::std::FUNC
#    define PONCA_MULTIARCH_INCLUDE_CU_STD(FILENAME) <cuda/std/FILENAME>

#    define PONCA_CUDA_ARCH
#else
#    define PONCA_MULTIARCH_INCLUDE_STD(FILENAME) <FILENAME>
#    define PONCA_MULTIARCH_CU_STD_FUNC(FUNC) using std::FUNC;
#    define PONCA_MULTIARCH_STD_MATH(FUNC) using std::FUNC;
#    define PONCA_MULTIARCH_STD_MATH_NAMESPACE(FUNC) std::FUNC
#    define PONCA_MULTIARCH_CU_STD_NAMESPACE(FUNC) std::FUNC
#    define PONCA_MULTIARCH_INCLUDE_CU_STD(FILENAME) <FILENAME>
#    define PONCA_CPU_ARCH
#endif

#ifndef PONCA_DEBUG
#    define STD_SAFE_AT(C, i) C[i]
#else
#    define STD_SAFE_AT(C, i) C.at(i)
#endif

#ifndef M_PI
// Source: http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
#    define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

namespace Ponca
{
    /*!
     * \brief Utility structure used to detect if a Point has a normal field
     *
     * Example usage:
     * \code
     * hasNormal<MyPoint>::value will be true if MyPoint has a member function 'normal()'
     * \endcode
     *
     * \tparam T The Point type
     */
    template <typename T, typename = void>
    struct hasNormal : std::false_type
    {
    };

    /// \copydoc hasNormal<typename,typename>
    template <typename T>
    struct hasNormal<T, std::void_t<decltype(std::declval<T>().normal())>> : std::true_type
    {
    };
} // namespace Ponca

