/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include <cstddef>
#include "../defines.h"
#include "./iteratorUtils.h"
#include "../../Common/Assert.h"

namespace Ponca
{
    /*! \brief A simple BitSet implementation that mimics a set of indices.
     * The internal logic of this bitset is similar to `std::bitset`, but is compatible with CUDA.
     *
     * Allows to insert and search in O(1) complexity, but is memory expensive, because we allocate a single
     * bit for each possible indices, to flag if it was inserted or not, which is not ideal for large amount of indices.
     *
     * The memory use of the BitSet depending on N should be in theory :
     *
     * BitSet size   | Memory used
     * ------------- | ------------
     * Bitset<10000> | 1.25 KB
     * Bitset<1e6>   | 125 KB
     * Bitset<1e7>   | 1.25 MB
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
        PONCA_MULTIARCH BitSet() = default;

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
        PONCA_MULTIARCH [[nodiscard]] bool contains(int value) const;

        //! \brief Sets all the bits to EMPTY
        PONCA_MULTIARCH void clear();

    protected:
        static constexpr size_t BIT_SIZE   = sizeof(T) * 8; //! The number of bits in one element of the array
        static constexpr size_t ARRAY_SIZE = (N + BIT_SIZE - 1) / BIT_SIZE; //!< The array size
        T m_data[ARRAY_SIZE]               = {};                            //!< An array of bytes
    };
} // namespace Ponca

#include "bitset.hpp"
