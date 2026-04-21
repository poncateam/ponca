#pragma once

#include <cstdio>

// macro delimiters
#define PONCA_MACRO_START \
    do                    \
    {
#define PONCA_MACRO_END \
    }                   \
    while (0)

// Stringification
#define PONCA_XSTR(S) #S
#define PONCA_STR(S) PONCA_XSTR(S)

#ifdef __CUDA_ARCH__
#    define PONCA_ABORT asm("trap;");
#elif defined(_MSC_VER)
#    define PONCA_ABORT __debugbreak();
#elif defined(__clang__) || defined(__GNUC__)
#    define PONCA_ABORT __builtin_trap();
#else
#    include <cstdlib>
#    define PONCA_ABORT std::abort();
#endif

#define PONCA_CRASH   \
    PONCA_MACRO_START \
    PONCA_ABORT       \
    PONCA_MACRO_END

#ifdef __CUDA_ARCH__
#    define PONCA_PRINT_ERROR(MSG)                              \
        PONCA_MACRO_START                                       \
        printf("%s:%i: [Error] %s\n", __FILE__, __LINE__, MSG); \
        PONCA_MACRO_END
#    define PONCA_PRINT_WARNING(MSG)                              \
        PONCA_MACRO_START                                         \
        printf("%s:%i: [Warning] %s\n", __FILE__, __LINE__, MSG); \
        PONCA_MACRO_END
#else
#    define PONCA_PRINT_ERROR(MSG)                                       \
        PONCA_MACRO_START                                                \
        fprintf(stderr, "%s:%i: [Error] %s\n", __FILE__, __LINE__, MSG); \
        fflush(stderr);                                                  \
        PONCA_MACRO_END
#    define PONCA_PRINT_WARNING(MSG)                                       \
        PONCA_MACRO_START                                                  \
        fprintf(stderr, "%s:%i: [Warning] %s\n", __FILE__, __LINE__, MSG); \
        fflush(stderr);                                                    \
        PONCA_MACRO_END
#endif

// turnoff warning
#define PONCA_UNUSED(VAR)         \
    PONCA_MACRO_START(void)(VAR); \
    PONCA_MACRO_END

#define PONCA_TODO             \
    PONCA_MACRO_START          \
    PONCA_PRINT_ERROR("TODO"); \
    PONCA_CRASH;               \
    PONCA_MACRO_END

#ifdef __has_builtin
#    if __has_builtin(__builtin_clz)
#        define PONCA_HAS_BUILTIN_CLZ 1
#    endif
#endif

#ifndef PONCA_HAS_BUILTIN_CLZ
#    define PONCA_HAS_BUILTIN_CLZ 0
#endif
