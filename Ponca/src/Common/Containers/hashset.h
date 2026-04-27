/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include "../defines.h"
#include "./iteratorUtils.h"

namespace Ponca
{
    /*! \brief A simple HashMap implementation that mimics a set of indices, by only storing the keys (Hence the name
     * HashSet). The internal logic of this HashMap is similar to a `std::unordered_map`, but is compatible with CUDA.
     *
     * Stores unique signed integer values in a contiguous array.
     *
     * \warning The internal values stored must not be equal to -1, because it is used as an empty flag by the HashSet.
     *
     * For the search, the best case complexity is O(1), and the worst case complexity is O(n).
     * The complexity is entirely dependent on the hashing function, and the given values :
     * If our hashing function produce sparser result for our set of indices, it will reduce the complexity of the
     * searches, and therefore improve the performance of our algorithm.
     *
     * \see HashSet::hash for the hashing function
     *
     * \tparam N The maximum size of the HashSet
     * \tparam T The value type stored in the HashSet : Must be a signed integer type, because the values are used as
     * indices for the array
     */
    template <int N, typename T = int>
    class HashSet
    {
        static_assert(N > 0, "The capacity must be strictly positive");

    public:
        /*! \brief Empty the array
         *
         * Iterates over every element and sets their values to -1 (empty flag)
         */
        PONCA_MULTIARCH void clear();

        /*! \brief Tries to insert a value in the HashSet
         *
         * \note The hashing function doesn't guarantee that we always get different results, sometimes two different
         * values can happen to have the same storage address in our data array. When the hashing function provokes a
         * collision with another value, the HashSet will store the other value in the next available array address. The
         * search will have to keep looking for the value if it isn't stored at the hashing result, which is why the
         * search isn't always of a O(1) complexity.
         *
         * \param value The value to be inserted in the HashSet
         * \return True if the value was inserted successfully, and false if the value was already inserted or if the
         * HashSet is full
         */
        PONCA_MULTIARCH bool insert(int value);

        /*! \brief Tries to find a value in the HashSet
         *
         * \param value The value to search for
         * \return True if the value is inside the HashSet, false if it's not in the HashSet.
         */
        PONCA_MULTIARCH bool contains(int value);

    private:
        static constexpr T EMPTY = T(-1);
        T m_data[N]              = { EMPTY };

        //! \brief The hashing function : (x * 2654435761u) % N
        PONCA_MULTIARCH [[nodiscard]] static int hash(const int x)
        {
            PONCA_ASSERT_MSG(x >= 0, "Index must be positive HashSet");
            return (x * 2654435761u) % N;
        }
    };

    template <int N, typename T>
    PONCA_MULTIARCH void HashSet<N, T>::clear()
    {
        Ponca::internal::fill(m_data, m_data + N, EMPTY);
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::insert(const int value)
    {
        const int h = hash(value);

        // Try to insert
        for (int i = 0; i < N; ++i)
        {
            const int idx = (h + i) % N;
            T& slot       = m_data[idx]; // Get the address

            // Stores here if the address is empty
            if (slot == EMPTY)
            {
                slot = value;
                return true;
            }
            // The value was already inserted in the array
            if (slot == value)
            {
                return false;
            }
        }

        // The array is full
        return false;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::contains(const int value)
    {
        const int h = hash(value);

        // Try to find the value
        for (int i = 0; i < N; ++i)
        {
            const int idx = (h + i) % N;
            T& slot       = m_data[idx]; // Get the address

            // Stops the search here if the address is empty
            if (slot == EMPTY)
            {
                return false;
            }
            // Value found
            if (slot == value)
            {
                return true;
            }
            // The value might have been inserted elsewhere, keep looking...
        }
        // Value not in HashSet
        return false;
    }
} // namespace Ponca
