/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <type_traits>
#include <string>

/**
 * \brief Mangles a Scalar
 *
 * \tparam T Scalar type
 */
template <typename T>
inline std::string MangleType()
{
    if constexpr (std::is_same_v<T, float>)
        return "f";
    if constexpr (std::is_same_v<T, double>)
        return "d";
    return "unknown";
}

/**
 * \brief Mangles an nb::ndarray::dtype
 *
 * This function is templated so that no constraints is
 * imposed on the actual nb::ndarray type.
 *
 * \tparam ArrType a nb::ndarray instance.
 */
template <typename ArrType>
inline std::string MangleDType(const ArrType& arr)
{
    namespace nb = nanobind;

    if (arr.dtype() == nb::dtype<double>())
        return MangleType<double>();
    if (arr.dtype() == nb::dtype<float>())
        return MangleType<float>();
    return "unknown";
}

/**
 * \brief Mangles an nb::ndarray
 *
 * This function is templated so that no constraints is
 * imposed on the actual nb::ndarray type.
 *
 * \tparam ArrType a nb::ndarray instance.
 */
template <typename ArrType>
inline std::string MangleArray(const ArrType& arr)
{
    std::string name = "";
    if (arr.ndim() >= 2)
        name = std::to_string(arr.shape(1));

    for (size_t i = 2; i < arr.ndim(); ++i)
        name += "_" + std::to_string(arr.shape(i));

    return name + MangleDType(arr);
}

/**
 * \brief Mangles a Ponca Point
 *
 * \tparam P Point type
 */
template <typename P>
inline std::string ManglePoint()
{
    return std::to_string(P::Dim) + MangleType<typename P::Scalar>();
}

