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
    template <int N, typename T = int>
    struct HashDefaultFunctor
    {
        //! \brief The default hashing function : (abs(x) * 2654435761u) % N
        PONCA_MULTIARCH [[nodiscard]] static constexpr T hash(const int _x)
        {
            PONCA_MULTIARCH_STD_MATH(abs);
            return (abs(_x) * 2654435761u) % N;
        }
    };

    /*! \brief A simple HashMap implementation that mimics a set of indices, by only storing the keys (Hence the name
     * HashSet). The internal logic of this HashMap is similar to a `std::unordered_map`, but is compatible with CUDA.
     *
     * Stores unique integer values in a contiguous array.
     *
     * \warning Logic was optimized for signed integer: The stored values must never be equal to -1, because it
     * will be mistaken as being empty in the HashSet, and break the search logic. Change the OFFSET value, depending
     * on the negative value you need to store to avoid this issue, with the following rule : illegal_value =
     * EMPTY-OFFSET (e.g. set OFFSET to 2 to allow to store -1, but make -2 illegal to store).
     *
     * For the search, the best case complexity is O(1), and the worst case complexity is O(n).
     * The complexity is entirely dependent on the hashing function, and the given values :
     * If our hashing function produce sparser result for our set of indices, it will reduce the complexity of the
     * searches, and therefore improve the performance of our algorithm.
     *
     * \see HashSet::hash for the hashing function
     * \see BitSet For alternative data structure with compatible API
     *
     * \tparam N The maximum size of the HashSet
     * \tparam T The value type stored in the HashSet. Default to int
     */
    template <int N, typename T = int, template <int, typename> typename _HashFunctor = HashDefaultFunctor>
    class HashSet
    {
        static_assert(N > 0, "The capacity must be strictly positive");
        using HashFunctor = _HashFunctor<N, T>;

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
         * \param _hash The hashing function
         * \return True if the value is inside the HashSet, false if it's not in the HashSet.
         */
        PONCA_MULTIARCH [[nodiscard]] inline bool search(int _value, int& _searchedIdx) const;

    public:
        constexpr PONCA_MULTIARCH HashSet() : m_data()
        {
            // Skip this initialization step if EMPTY is set to 0
            if constexpr (EMPTY != T(0))
            {
                Ponca::internal::fill(m_data, m_data + N, EMPTY);
            }
        }

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
         * \param _value The value to be inserted in the HashSet
         * \return True if the value was inserted successfully, and false if the value was already inserted or if the
         * HashSet is full
         */
        PONCA_MULTIARCH bool insert(int _value);

        /*! \brief Tries to find a value in the HashSet
         *
         * \param _value The value to search for
         * \return True if the value is inside the HashSet, false if it's not in the HashSet.
         */
        PONCA_MULTIARCH [[nodiscard]] bool contains(int _value) const;

    private:
        static constexpr T OFFSET =
            T(1); //< Offsets the value when storing in m_data, to avoid mistaking the stored index value with EMPTY
        static constexpr T EMPTY = T(0); //< The flag to tell if the address is available or not (Should always be zero)
        T m_data[N];                     //< Where we store the elements in memory
    };
} // namespace Ponca

#include "./hashset.hpp"
