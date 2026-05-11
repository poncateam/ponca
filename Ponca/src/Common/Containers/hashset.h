/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include "../defines.h"
#include "./iteratorUtils.h"
#include <utility>

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

    /*! \brief Stores unique signed integer values in a contiguous array.
     *
     * A simple HashMap implementation that mimics a set of indices, by only stores the keys in a Set-like
     * structure (Hence the name HashSet). The internal logic of this HashMap is similar to a `std::unordered_map`,
     * but is compatible with CUDA.
     *
     * \note This HashSet type doesn't allow for removal of a single element from the set. It was implemented this way
     * to provide the simplest search and insert method possible
     *
     * \warning The stored values must never be equal to -1, because it will be mistaken as being empty in the HashSet,
     * and will break the search logic. If you must store -1, change the OFFSET value to something else, by following
     * this simple rule : illegal_value = -OFFSET (e.g. set OFFSET to 2 to allow to store -1, but make -2 illegal
     * to store).
     *
     * For the search, the best case complexity is O(1), and the worst case complexity is O(n).
     * The complexity is entirely dependent on the hashing function, and the given values :
     * If our hashing function produce sparser result for our set of indices, it will reduce the complexity of the
     * searches for both the `HashSet::insert` and `HashSet::contains` method.
     *
     * \see HashDefaultFunctor::hash for the default hashing function
     * \see BitSet For alternative data structure with compatible API
     *
     * \tparam N The maximum size of the HashSet
     * \tparam T The value type stored in the HashSet (Default to int)
     * \tparam _HashFunctor For the hashing function (Default to HashDefaultFunctor)
     * \tparam OFFSET Offsets the value stored in m_data, to avoid mistaking it with 0 the is empty flag
     */
    template <int N, typename T = int, template <int, typename> typename _HashFunctor = HashDefaultFunctor,
              T OFFSET = T(1)>
    class HashSet
    {
        static_assert(N > 0, "The capacity must be strictly positive");
        using HashFunctor    = _HashFunctor<N, T>;
        using container_type = std::array<T, N>;
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        using Self           = HashSet<N, T, _HashFunctor, OFFSET>;

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
        PONCA_MULTIARCH [[nodiscard]] inline bool search(T _value, T& _searchedIdx) const;

    public:
        constexpr PONCA_MULTIARCH HashSet() : m_data() {}

        /*! \brief Empty the array
         *
         * Iterates over every element and sets their values to 0
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
         * \note Throws an error if value was inserted in an already full HashSet
         *
         * \param _value The value to be inserted in the HashSet
         * \return An std::pair of
         * - First : The iterator to the inserted value or end() if the HashSet is full
         * - Second : True if the value was inserted successfully, and false if the value wasn't inserted
         */
        PONCA_MULTIARCH std::pair<typename Self::iterator, bool> insert(const T& _value);

        /*! \brief Tries to find a value in the HashSet
         *
         * \param _value The value to search for
         * \return True if the value is inside the HashSet, false if it's not in the HashSet.
         */
        PONCA_MULTIARCH [[nodiscard]] bool contains(T _value) const;

    public:
        //! \brief The beginning of the internal array
        PONCA_MULTIARCH [[nodiscard]] inline Self::const_iterator cbegin() const;

        //! \brief The end of the internal array
        PONCA_MULTIARCH [[nodiscard]] inline Self::const_iterator cend() const;

        //! \brief The beginning of the internal array
        PONCA_MULTIARCH [[nodiscard]] inline Self::iterator begin();

        //! \brief The end of the internal array
        PONCA_MULTIARCH [[nodiscard]] inline Self::iterator end();

    private:
        container_type m_data{}; //< Where we store the elements in memory
    };
} // namespace Ponca

#include "./hashset.hpp"
