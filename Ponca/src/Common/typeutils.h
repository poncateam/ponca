/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

namespace Ponca
{
    /**
     * \brief StringLiteral as template parameters
     *
     * We use the C++20 NTTP extension to allow to template by string literals.
     * This class will be used to provide a name to FactoryEntries
     */
    template <size_t N>
    struct StringLiteral
    {
        char value[N];
        constexpr StringLiteral(const char (&str)[N]) { std::copy_n(str, N, value); }
    };

    /**
     * \brief Constexpr version of string comparison
     *
     * At the time of writting this function, strcmp is not constexpr
     * and hence can not be used in templated code.
     *
     * \param a First string
     * \param b Second string
     */
    inline constexpr bool strequal(const char* a, const char* b)
    {
        while (*a && *b)
        {
            if (*a != *b)
                return false;

            ++a;
            ++b;
        }

        return *a == *b;
    }
} // namespace Ponca
