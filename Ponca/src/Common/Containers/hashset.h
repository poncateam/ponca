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

    protected:
        /*! \brief Search for a value in the HashSet.
         *
         * Use the hashing function to check if a value is inside the internal array.
         *
         * - The value isn't considered inside if the value at the address corresponds to the EMPTY flag, (will stop
         * the search and return false as a result)
         *
         * - If there is already an element at address, but which doesn't correspond to the value, we keep looking at
         * the next address for our value until we either find it, or until we searched everywhere.
         *
         * Use a reference to either output the id of the last element that was searched in the m_data array or will
         * output -1 if the array is completely full.
         *
         * \param _value The value to be inserted in the HashSet
         * \param _searchedIdx Reference to the last searched index or -1 if the array is full.
         * \return True if the value is inside the HashSet, false if it's not in the HashSet.
         */
        PONCA_MULTIARCH [[nodiscard]] bool search(int _value, int& _searchedIdx) const;

    public:
        /*! \brief Empty the array
         *
         * Iterates over every element and sets their values to the EMPTY flag value (default to -1)
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
        PONCA_MULTIARCH [[nodiscard]] bool contains(int value) const;

    private:
        static constexpr T EMPTY = T(-1);
        T m_data[N]              = {EMPTY};

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
    PONCA_MULTIARCH bool HashSet<N, T>::search(const int _value, int& _searchedIdx) const
    {
        const int h = hash(_value);

        // Try to find the value
        for (int i = 0; i < N; ++i)
        {
            _searchedIdx  = (h + i) % N;
            const T& slot = m_data[_searchedIdx]; // Get the address

            // Stops the search here if the address is empty
            if (slot == EMPTY)
                return false;

            // Value found
            if (slot == _value)
                return true;

            // The value might have been inserted elsewhere, keep looking...
        }

        // Value not in HashSet
        _searchedIdx = -1;
        return false;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::insert(const int value)
    {
        int availableIdx = 0;
        if (search(value, availableIdx))
        { // If search is successful
            // Insertion can't be done because found the value in the array
            std::cout << "Insertion can't be done because found duplicated at" << availableIdx << "." << std::endl;
            return false;
        }
        // The value wasn't found in the array, so either
        // A - The array is full (The last search index shouldn't point to an available address in the array)
        if (availableIdx == -1) // Search returns -1 if it was full
        {
            std::cout << "Insertion can't be done because array is full" << std::endl;
            return false;
        }

        // B - The array isn't full and the value can be inserted (The last search index should therefore point to an
        // available address in the array)
        m_data[availableIdx] = value;
        return true;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::contains(const int value) const
    {
        int i = 0;
        return search(value, i);
    }
} // namespace Ponca
