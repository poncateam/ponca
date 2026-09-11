/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <type_traits>

#include <Ponca/src/Common/typeutils.h>

/**
 * \brief Class to link a neighbor filter and its mangled name. 
 * 
 * We store the name into a template parameter because the current
 * paradigm does not use isntance but rather template types. 
 */
template<typename _NF, Ponca::StringLiteral _Name>
struct MangledFilter
{
    using NF = _NF;
    static constexpr Ponca::StringLiteral _name = _Name;
    static constexpr const char* name = _name.value;
};

// We use accronyms here to avoid name collision
template <typename Point>
using SWFilter = MangledFilter<Ponca::DistWeightFilter<Point, Ponca::SmoothWeightKernel<typename Point::Scalar>>, "SW">;

template <typename Point>
using CWFilter = MangledFilter<Ponca::DistWeightFilter<Point, Ponca::ConstantWeightKernel<typename Point::Scalar>>, "CW">;

template <typename Point>
using NWFilter = MangledFilter<Ponca::NoWeightFilter<Point>, "NW">;