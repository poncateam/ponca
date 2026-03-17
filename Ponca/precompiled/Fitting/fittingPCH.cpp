#include "fittingETI.h"

#define EXTERN
#include "fittingDeclareMacro.h"

/*!
  \file fittingPCH.cpp
  \brief Precompiled header target that defines the commonly used Fitting types for 3D point clouds.
*/

// Explicit template instantiation
FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef FITTING_DEF
#undef EXTERN
