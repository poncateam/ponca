/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

namespace Ponca
{
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

            // Is stored as value+OFFSET in the array (see insert)
            if (slot == _value + OFFSET)
                return true; // Value was found

            // The value might have been inserted elsewhere, keep looking...
        }

        // Value not in HashSet
        _searchedIdx = -1;
        return false;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::insert(const int _value)
    {
        PONCA_ASSERT_MSG(_value != EMPTY + OFFSET, "Illegal value was inserted into the HashSet");
        int availableIdx = 0;
        if (search(_value, availableIdx)) // If search is successful
            return false;                 // Insertion can't be done because found the value in the array

        // The value wasn't found in the array, so either :
        // A - The set is full (The last search index shouldn't point to an available address in the array)
        if (availableIdx == -1) // Search returns -1 if the Set is full
            return false;

        // B - The set isn't full and the value can be inserted. Therefore, the last search index is the next available
        // address in the array
        m_data[availableIdx] = _value + OFFSET;
        return true;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::contains(const int _value) const
    {
        PONCA_DEBUG_ASSERT_MSG(_value != EMPTY + OFFSET, "Illegal value was searched from the HashSet");
        int i = 0;
        return search(_value, i);
    }
} // namespace Ponca
