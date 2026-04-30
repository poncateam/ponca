/**
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble - Main source code
 \author Nicolas Mellado   - Improved error handling and bound check
*/

namespace Ponca
{
    template <class T, int N>
    Stack<T, N>::Stack() : m_size(0), m_data()
    {
    }

    template <class T, int N>
    const T& Stack<T, N>::top() const
    {
        return STD_SAFE_AT(m_data, m_size - 1);
    }

    template <class T, int N>
    T& Stack<T, N>::top()
    {
        return STD_SAFE_AT(m_data, m_size - 1);
    }

    template <class T, int N>
    bool Stack<T, N>::empty() const
    {
        return m_size == 0;
    }

    template <class T, int N>
    int Stack<T, N>::size() const
    {
        return m_size;
    }

    template <class T, int N>
    void Stack<T, N>::push(const T& value)
    {
        STD_SAFE_AT(m_data, m_size) = value;
        ++m_size;
    }

    template <class T, int N>
    void Stack<T, N>::push()
    {
        ++m_size;
    }

    template <class T, int N>
    void Stack<T, N>::pop()
    {
        --m_size;
    }

    template <class T, int N>
    void Stack<T, N>::clear()
    {
        m_size = 0;
    }

} // namespace Ponca
