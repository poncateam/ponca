/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#ifdef _PONCA_COMPILE_DEFINITION
#    define _PRECOMPILED_PONCA_EXTERN
#else
#    define _PRECOMPILED_PONCA_EXTERN extern
#endif

#include "../../Ponca"

namespace Ponca
{
// Note: Variadic macros here because ',' in template definition are parsed as different arguments
// we could use a trick by wrapping the class in parenthesis, but for now this works and is simpler.
#define _PONCA_BASKET_X(name, desc, ...) _PRECOMPILED_PONCA_EXTERN template class __VA_ARGS__;
#define _PONCA_BASKET_DIFF_X(name, desc, ...) _PRECOMPILED_PONCA_EXTERN template class __VA_ARGS__;
#include "instantiate.h"
#undef _PONCA_BASKET_X
#undef _PONCA_BASKET_DIFF_X
}; // namespace Ponca
