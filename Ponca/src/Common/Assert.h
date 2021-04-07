#pragma once

#include "./Macro.h"

#define PONCA_ERROR                                                             \
    PONCA_MACRO_START                                                           \
    PONCA_PRINT_ERROR("");                                                      \
    PONCA_CRASH;                                                                \
    PONCA_MACRO_END

#define PONCA_ERROR_MSG(MSG)                                                    \
    PONCA_MACRO_START                                                           \
    PONCA_PRINT_ERROR(MSG);                                                     \
    PONCA_CRASH;                                                                \
    PONCA_MACRO_END

#define PONCA_ASSERT(EXPR)                                                      \
    PONCA_MACRO_START                                                           \
    if(!(EXPR)) {                                                              \
        PONCA_ERROR_MSG(PONCA_STR(EXPR));                                        \
    }                                                                          \
    PONCA_MACRO_END

#define PONCA_ASSERT_MSG(EXPR,MSG)                                              \
    PONCA_MACRO_START                                                           \
    if(!(EXPR)) {                                                              \
        PONCA_ERROR_MSG(MSG);                                                   \
    }                                                                          \
    PONCA_MACRO_END

#ifdef PONCA_DEBUG
    #define PONCA_DEBUG_ASSERT(EXPR)          PONCA_ASSERT(EXPR)
    #define PONCA_DEBUG_ASSERT_MSG(EXPR,MSG)  PONCA_ASSERT_MSG(EXPR,MSG)
    #define PONCA_DEBUG_ERROR_MSG(MSG)        PONCA_ERROR_MSG(MSG)
    #define PONCA_DEBUG_ERROR                 PONCA_ERROR
#else
    #define PONCA_DEBUG_ASSERT(EXPR)
    #define PONCA_DEBUG_ASSERT_MSG(EXPR,MSG)
    #define PONCA_DEBUG_ERROR_MSG(MSG)
    #define PONCA_DEBUG_ERROR
#endif
