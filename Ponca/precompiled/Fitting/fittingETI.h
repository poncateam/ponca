/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file Ponca/precompiled/Fitting/fittingETI.h
 * \brief Pre-instantiate explicitly most fitting types for 3D point clouds.
 * It tells the compiler that it should fetch with the vtable the Fitting types from the precompiled library.
 */

#pragma once

#define EXTERN extern
#include "fittingDeclareMacro.h"

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef EXTERN
