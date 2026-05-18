/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

namespace Ponca
{
    // Iterators --------------------------------------------------------------------
    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    typename HashSet<N, T, HF, OFFSET>::const_iterator HashSet<N, T, HF, OFFSET>::cbegin() const
    {
        return m_data.begin();
    }

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    typename HashSet<N, T, HF, OFFSET>::const_iterator HashSet<N, T, HF, OFFSET>::cend() const
    {
        return m_data.begin() + N;
    }

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    typename HashSet<N, T, HF, OFFSET>::iterator HashSet<N, T, HF, OFFSET>::begin()
    {
        return m_data.begin();
    }

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    typename HashSet<N, T, HF, OFFSET>::iterator HashSet<N, T, HF, OFFSET>::end()
    {
        return m_data.begin() + N;
    }

    // Set-like methods --------------------------------------------------------------------
    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    void HashSet<N, T, HF, OFFSET>::clear()
    {
        Ponca::internal::fill(m_data.begin(), m_data.begin() + N, 0);
    }

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    bool HashSet<N, T, HF, OFFSET>::search(const T _value, T& _searchedIdx) const
    {
        const int h = HashFunctor::hash(_value);

        // Try to find the value
        for (int i = 0; i < N; ++i)
        {
            _searchedIdx  = (h + i) % N;
            const T& slot = m_data[_searchedIdx]; // Get the address

            // Stops the search here if the address is empty
            if (slot == 0)
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

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    std::pair<typename HashSet<N, T, HF, OFFSET>::iterator, bool> HashSet<N, T, HF, OFFSET>::insert(const T& _value)
    {
        PONCA_ASSERT_MSG(_value != -OFFSET, "Illegal value was inserted into the HashSet");
        int availableIdx = 0;
        if (search(_value, availableIdx)) // If search is successful
            return std::make_pair(m_data.begin() + availableIdx,
                                  false); // Insertion can't be done because found the value in the array

        // The value wasn't found in the array, so either :
        // A - The set is full (The last search index is -1 if it didn't find an available address in the array)
        if (availableIdx == -1)
        {
            return std::make_pair(end(), false);
        }
        // B - The set isn't full and the value can be inserted. Therefore, the last search index is the next available
        // address in the array
        m_data[availableIdx] = _value + OFFSET;
        return std::make_pair(m_data.begin() + availableIdx, true);
    }

    template <int N, typename T, template <int, typename> typename HF, T OFFSET>
    bool HashSet<N, T, HF, OFFSET>::contains(T _value) const
    {
        PONCA_DEBUG_ASSERT_MSG(_value != -OFFSET, "Illegal value was searched from the HashSet");
        int i;
        return search(_value, i);
    }
} // namespace Ponca
