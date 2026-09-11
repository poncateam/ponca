/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace Ponca
{
    // TODO: Specialized bitwise version for 32/64 bits

    template <unsigned int D, typename Int = std::uint32_t>
    constexpr inline Int MortonEncode(const std::array<Int, D>& coords)
    {
        constexpr size_t bits = sizeof(Int) * 8 / D;

        Int code = 0;
        for (size_t b = 0; b < bits; ++b)
        {
            for (size_t d = 0; d < D; ++d)
            {
                // Get b-th bit and put it in position b * D + d
                code |= ((coords[d] >> b) & Int{1}) << (b * D + d);
            }
        }
        return code;
    }

    template <unsigned int D, typename Int = std::uint32_t>
    constexpr inline std::array<Int, D> MortonDecode(Int code)
    {
        constexpr size_t bits = sizeof(Int) * 8 / D;

        std::array<Int, D> coords{};
        for (size_t b = 0; b < bits; ++b)
        {
            for (size_t d = 0; d < D; ++d)
            {
                coords[d] |= ((code >> (b * D + d)) & Int{1}) << b;
            }
        }
        return coords;
    }
} // namespace Ponca
