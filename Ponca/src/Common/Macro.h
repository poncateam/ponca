#pragma once

#include <cstdio>

// macro delimiters
#define PONCA_MACRO_START do {
#define PONCA_MACRO_END   } while(0)

// Stringification
#define PONCA_XSTR(S) #S
#define PONCA_STR(S) PONCA_XSTR(S)

#define PONCA_CRASH                                                             \
    PONCA_MACRO_START                                                           \
    asm volatile ("int $3");                                                   \
    PONCA_MACRO_END

#define PONCA_PRINT_ERROR(MSG)                                                  \
    PONCA_MACRO_START                                                           \
    fprintf(stderr,                                                            \
            "%s:%i: [Error] %s\n",                                             \
            __FILE__,__LINE__,MSG);                                            \
    fflush(stderr);                                                            \
    PONCA_MACRO_END

#define PONCA_PRINT_WARNING(MSG)                                                \
    PONCA_MACRO_START                                                           \
    fprintf(stderr,                                                            \
            "%s:%i: [Warning] %s\n",                                           \
            __FILE__,__LINE__,MSG);                                            \
    fflush(stderr);                                                            \
    PONCA_MACRO_END

// turnoff warning
#define PONCA_UNUSED(VAR)                                                       \
    PONCA_MACRO_START                                                           \
    (void)(VAR);                                                               \
    PONCA_MACRO_END

#define PONCA_TODO                                                              \
    PONCA_MACRO_START                                                           \
    PONCA_PRINT_ERROR("TODO");                                                  \
    PONCA_CRASH;                                                                \
    PONCA_MACRO_END

template<typename T> struct XXXXXXXXXX_The_unkown_type_is_;
#define WhatIsTheTypeOf(expr) XXXXXXXXXX_The_unkown_type_is_<decltype(expr)> _XXXXXXXXXX
#define WhatIsThisType(type) XXXXXXXXXX_The_unkown_type_is_<type> _XXXXXXXXXX
