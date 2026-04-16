/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include "../defines.h"

namespace Ponca {
    /*! \brief A simple HashMap implementation, that stores unique elements in an array of signed values.
     *
     * The internal values stored must not be equal to -1, because it is used as an empty flag by the HashSet.
     * For the search, the best case complexity is O(1), and the worst case complexity is O(n).
     * The complexity is entirely dependent on the hashing function, and the given values :
     * Sparser hashing results will reduce the searching complexity.
     *
     * \tparam N The maximum size of the HashSet
     * \tparam T The value type stored in the HashSet : Must be a signed value type that can be used as an array index
     * (int-like)
     */
    template <int N, typename T>
    class HashSet {
    public:
        /*! \brief Empty the array
         *
         * Iterates over every element and sets its value to -1
         */
        PONCA_MULTIARCH void clear();

        /*! \brief Tries to insert a value in the HashSet
         *
         * When the hashing function provokes a collision (if two values have the same storage address in the data
         * array), the HashSet will store the value in the next available array address.
         * \param value The value to be inserted in the HashSet
         * \return True if the value was inserted successfully, and false if the value was already inserted or if the
         * HashSet is full
         */
        PONCA_MULTIARCH bool insert(int value);

    private:
        static constexpr T EMPTY = T(-1);
        T table[N] = {};

        PONCA_MULTIARCH [[nodiscard]] static int hash(const int x)
        {
            return (x * 2654435761u) % N;
        }
    };

    template <int N, typename T>
    PONCA_MULTIARCH void HashSet<N, T>::clear()
    {
        for (int i = 0; i < N; i ++)
            table[i] = EMPTY;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::insert(const int value)
    {
        const int h = hash(value);

        // Try to insert
        for (int i = 0; i < N; ++i) {
            int idx = (h + i) % N;

            int& slot = table[idx];

            // Stores here if the address is empty
            if (slot == EMPTY) {
                slot = value;
                return true;
            }
            // The value was already inserted in the array
            if (slot == value) {
                return false;
            }
        }

        // The array is full
        return false;
    }
}