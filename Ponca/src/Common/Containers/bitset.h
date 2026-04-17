/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include <cstddef>
#include "../defines.h"
#include <Ponca/src/Common/Assert.h>

namespace Ponca
{
    /*! \brief BitSet implementation, similar to std::bitset
     *
     * Allows to insert and search in O(1) complexity, but is memory expensive, because we allocate a single
     * bit for each possible indices, to flag if it was inserted or not, which is not ideal for large amount of indices.
     *
     * The memory use should be :
     * Bitset<10000> → 10k bits = 1.25 KB
     * Bitset<1e6> → 125 KB
     * Bitset<1e7> → 1.25 MB
     *
     * \tparam N Maximum number of indices
     * \tparam T The data type of the array storing the bits. Default to unsigned long long for 64 bits storage.
     */
    template <int N, typename T=unsigned long long>
    class BitSet
    {
        static_assert(N > 0, "The capacity must be strictly positive");
    public:
        PONCA_MULTIARCH BitSet();
        PONCA_MULTIARCH void flip(int index);
        PONCA_MULTIARCH bool insert(int index);
        PONCA_MULTIARCH [[nodiscard]] bool find(int index) const;
        PONCA_MULTIARCH void clear();

    protected:
        static constexpr size_t BIT_SIZE = sizeof(T) * 8; //! The number of bits in one element of the array
        static constexpr size_t ARRAY_SIZE = (N + BIT_SIZE - 1) / BIT_SIZE; //!< The array size
        T m_data[ARRAY_SIZE] = {}; //!< An array of bytes
    };

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    template <int N, typename T>
    BitSet<N, T>::BitSet() = default;

    template <int N, typename T>
    PONCA_MULTIARCH void BitSet<N, T>::clear()
    {
        for (int i = 0; i < ARRAY_SIZE; i++)
        {
            m_data[i] = T(0);
        }
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool BitSet<N, T>::insert(const int index) {
        PONCA_DEBUG_ASSERT(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        const T bitMask = (T(1) << bit);
        const bool alreadyInserted = (m_data[byte] & bitMask) == 0;
        m_data[byte] |= bitMask;
        return alreadyInserted;
    }

    template <int N, typename T>
    PONCA_MULTIARCH void BitSet<N, T>::flip(const int index) {
        PONCA_DEBUG_ASSERT(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        m_data[byte] ^= (T(1) << bit);
    }

    template <int N, typename T>
    PONCA_MULTIARCH [[nodiscard]] bool BitSet<N, T>::find(const int index) const {
        PONCA_DEBUG_ASSERT(index>=0 && index<N);
        const int byte = index / BIT_SIZE;
        const int bit  = index % BIT_SIZE;
        return (m_data[byte] & (T(1) << bit)) != 0;
    }
}
