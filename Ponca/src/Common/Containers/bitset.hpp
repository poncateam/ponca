/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

namespace Ponca
{
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// Set like methods ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    template <int N, typename T>
    PONCA_MULTIARCH void BitSet<N, T>::clear()
    {
        Ponca::internal::fill(m_data, m_data + ARRAY_SIZE, T(0));
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
    PONCA_MULTIARCH [[nodiscard]] bool BitSet<N, T>::contains(const int value) const
    {
        PONCA_ASSERT_MSG(value >= 0 && value < N, "Searched value is outside the scope of the BitSet");
        const int byte = value / BIT_SIZE;
        const int bit  = value % BIT_SIZE;
        return (m_data[byte] & (T(1) << bit)) != 0;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// BitSet like methods ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    template <int N, typename T>
    PONCA_MULTIARCH void BitSet<N, T>::flip(const int i)
    {
        PONCA_ASSERT_MSG(i >= 0 && i < N, "Flipped value is outside the scope of the BitSet");
        const int byte = i / BIT_SIZE;
        const int bit  = i % BIT_SIZE;
        m_data[byte] ^= (T(1) << bit);
    }
} // namespace Ponca