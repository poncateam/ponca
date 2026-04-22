/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include <cstddef>
#include "../defines.h"

namespace Ponca
{
    /*! \brief A simple BitSet implementation that mimics a set of indices.
     * The internal logic of this bitset is similar to `std::bitset`, but is compatible with CUDA.
     *
     * Allows to insert and search in O(1) complexity, but is memory expensive, because we allocate a single
     * bit for each possible indices, to flag if it was inserted or not, which is not ideal for large amount of indices.
     *
     * The memory use of the BitSet depending on N should be :
     * Bitset<10000> → 10k bits = 1.25 KB
     * Bitset<1e6>   → 125 KB
     * Bitset<1e7>   → 1.25 MB
     *
     * \tparam N Maximum number of indices
     * \tparam T The data type of the array storing the bits. Default to 'unsigned long long' for 64 bits storage.
     *
     * \see `BitSet::erase`, `BitSet::insert` for set-index-like methods
     *
     * \warning The inserted values must always be smaller than N, because we are using them as indices to store the
     * bits inside the BitSet.
     */
    template <int N, typename T = unsigned long long>
    class BitSet
    {
        static_assert(N > 0, "The capacity must be strictly positive");

    public:
        PONCA_MULTIARCH BitSet();

        /*! \brief Toggles the value of a bit
         * \param i The bit to flip
         */
        PONCA_MULTIARCH void flip(int i);

        /*! \brief Tries to insert a value in the set
         *
         * \param value The value to be removed in the set. Must always be smaller than N, because we are using the
         * values as indices inside the BitSet.
         * \return True if the value was removed successfully, and false if the value was already not inside the Set
         */
        PONCA_MULTIARCH bool erase(int value);

        /*! \brief Tries to insert a value in the set
         *
         * \param value The value to be inserted in the HashSet. Must always be smaller than N, because we are using the
         * values as indices inside the BitSet.
         * \return True if the value was inserted successfully, and false if the value was already inserted or if the
         * HashSet is full
         */
        PONCA_MULTIARCH bool insert(int value);

        /*! \brief Search if the value was already inserted or not
         *
         * \param value The value to search for (search is O(N) because it corresponds to the index in the BitSet)
         * \return
         */
        PONCA_MULTIARCH [[nodiscard]] bool find(int value) const;

        //! \brief Sets all the bits to EMPTY
        PONCA_MULTIARCH void clear();

    protected:
        static constexpr size_t BIT_SIZE   = sizeof(T) * 8; //! The number of bits in one element of the array
        static constexpr size_t ARRAY_SIZE = (N + BIT_SIZE - 1) / BIT_SIZE; //!< The array size
        T m_data[ARRAY_SIZE]               = {};                            //!< An array of bytes
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
    PONCA_MULTIARCH bool BitSet<N, T>::erase(const int value)
    {
        PONCA_ASSERT_MSG(value >= 0 && value < N,
                         "Attempted to remove a value that is outside the scope of the BitSet");
        const int byte          = value / BIT_SIZE;
        const int bit           = value % BIT_SIZE;
        const T bitMask         = (T(1) << bit);
        const bool alreadyEmpty = (m_data[byte] & bitMask) != 0;
        m_data[byte] &= ~bitMask;
        return alreadyEmpty;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool BitSet<N, T>::insert(const int value)
    {
        PONCA_ASSERT_MSG(value >= 0 && value < N, "Inserted value is outside the scope of the BitSet");
        const int byte             = value / BIT_SIZE;
        const int bit              = value % BIT_SIZE;
        const T bitMask            = (T(1) << bit);
        const bool alreadyInserted = (m_data[byte] & bitMask) == 0;
        m_data[byte] |= bitMask;
        return alreadyInserted;
    }

    template <int N, typename T>
    PONCA_MULTIARCH void BitSet<N, T>::flip(const int i)
    {
        PONCA_ASSERT_MSG(i >= 0 && i < N, "Flipped value is outside the scope of the BitSet");
        const int byte = i / BIT_SIZE;
        const int bit  = i % BIT_SIZE;
        m_data[byte] ^= (T(1) << bit);
    }

    template <int N, typename T>
    PONCA_MULTIARCH [[nodiscard]] bool BitSet<N, T>::find(const int value) const
    {
        PONCA_ASSERT_MSG(value >= 0 && value < N, "Searched value is outside the scope of the BitSet");
        const int byte = value / BIT_SIZE;
        const int bit  = value % BIT_SIZE;
        return (m_data[byte] & (T(1) << bit)) != 0;
    }
} // namespace Ponca
