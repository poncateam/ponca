/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>

/**
 * \brief Register class, functions and enums for the common module
 * 
 * \param m The main module
 * \param internal An internal module reserved for the library
 */
void RegisterCommon(nanobind::module_& m, nanobind::module_& internal);
