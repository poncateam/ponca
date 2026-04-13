/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file Ponca/precompiled/Fitting/fittingPCH.h
 * \brief Link most fitting types to pre-instantiated class for 3D point clouds.
 * Tells the compiler that it should fetch the class using the vtable from the precompiled library.
 * This header is to be included in each target that wants to benefit from the pre-instantiation acceleration.
 */

#pragma once

#define EXTERN extern // Link to the vtable with the extern keyword
#include "fittingDeclareMacro.h"

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef EXTERN
