#pragma once

#define EXTERN extern
#include "fittingDeclareMacro.h"

/*!
  \file fittingETI.h
  \brief Predeclare the Explicit Template Instantiations that exist in the precompiled library.
  This header is to be included in each target that wants to benefit from precompiled Fitting acceleration.
  It tells the compiler that it should fetch the Fitting types from the precompiled library.
*/

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef EXTERN
