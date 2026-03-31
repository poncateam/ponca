/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include <concepts>

namespace Ponca
{
    template <typename T>
    concept ProvidesCommonTypes = requires(T t) {
        typename T::Scalar;
        typename T::VectorType;
    };

    template <typename T>
    concept IsPoint = ProvidesCommonTypes<T> && requires(const T ct) { ct.pos(); };

    template <typename T>
    concept IsPointNormal = IsPoint<T> && requires(const T ct) { ct.normal(); };

}; // namespace Ponca
