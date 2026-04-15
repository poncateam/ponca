/**
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include "../defines.h"

namespace Ponca
{
    //!
    //! @tparam N
    template <int N, typename T=unsigned long long>
    class Bitset
    {
    public:
        PONCA_MULTIARCH Bitset();
        PONCA_MULTIARCH void flip(int index);
        PONCA_MULTIARCH bool insert(int index);
        PONCA_MULTIARCH [[nodiscard]] bool find(int index) const;
        PONCA_MULTIARCH void clear();

    protected:
        static constexpr size_t BIT_SIZE = sizeof(T) * 8;
        static constexpr size_t ARRAY_SIZE = (N + BIT_SIZE - 1) / BIT_SIZE;
        T m_data[ARRAY_SIZE]; //!< An array of bytes
    };

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    template <int N, typename T>
    Bitset<N, T>::Bitset() = default;

    template <int N, typename T>
    PONCA_MULTIARCH void Bitset<N, T>::clear()
    {
        for (int i = 0; i < ARRAY_SIZE; i++)
        {
            m_data[i] = T(0);
        }
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool Bitset<N, T>::insert(const int index) {
        assert(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        const bool res = (m_data[byte] & (T(1) << bit)) != 0;
        m_data[byte] |= (T(1) << bit);
        return res;
    }

    template <int N, typename T>
    PONCA_MULTIARCH void Bitset<N, T>::flip(const int index) {
        assert(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        m_data[byte] ^= (T(1) << bit);
    }

    template <int N, typename T>
    PONCA_MULTIARCH [[nodiscard]] bool Bitset<N, T>::find(const int index) const {
        assert(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        return (m_data[byte] & (T(1) << bit)) != 0;
    }
}
