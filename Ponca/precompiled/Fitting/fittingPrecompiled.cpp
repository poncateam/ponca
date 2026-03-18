#include "fittingPCH.h"

/*!
  \file fittingPrecompiled.cpp
  \brief Precompiled target that defines the commonly used Fitting types for 3D point clouds.
*/

// Explicit template instantiation
#define EXTERN
#include "fittingDeclareMacro.h"

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef FITTING_DEF
#undef EXTERN
