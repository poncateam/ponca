/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

/**
  *
  * \defgroup fitting Fitting module
  * \brief This modules includes classes and methods for primitive fitting. See reference manual below.
  * @{
  * @}
  *
  */
#include "../Common/defines.h"

#define PONCA_EXPLICIT_CAST_OPERATORS(CLASSNAME,CONVERTER)                                                         \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
CLASSNAME<DataPoint, _WFunctor, T>& CONVERTER()                                                                    \
{ return * static_cast<CLASSNAME<DataPoint, _WFunctor, T>*>(this); }                                               \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
const CLASSNAME<DataPoint, _WFunctor, T>& CONVERTER() const                                                        \
{ return * static_cast<const CLASSNAME<DataPoint, _WFunctor, T>*>(this); }

#define PONCA_EXPLICIT_CAST_OPERATORS_DER(CLASSNAME,CONVERTER)                                                     \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
CLASSNAME<DataPoint, _WFunctor, DiffType, T>& CONVERTER()                                                          \
{ return * static_cast<CLASSNAME<DataPoint, _WFunctor, DiffType, T>*>(this); }                                     \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
const CLASSNAME<DataPoint, _WFunctor, DiffType, T>& CONVERTER() const                                              \
{ return * static_cast<const CLASSNAME<DataPoint, _WFunctor, DiffType, T>*>(this); }
