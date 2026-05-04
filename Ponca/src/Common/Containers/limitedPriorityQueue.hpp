/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble
*/

namespace Ponca
{
    // Iterator --------------------------------------------------------------------

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::iterator LimitedPriorityQueue<T, N, Cmp>::begin()
    {
        return m_data.begin();
    }

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::begin() const
    {
        return m_data.begin();
    }

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::cbegin() const
    {
        return m_data.cbegin();
    }

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::iterator LimitedPriorityQueue<T, N, Cmp>::end()
    {
        return m_data.begin() + m_size;
    }

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::end() const
    {
        return m_data.begin() + m_size;
    }

    template <class T, int N, class Cmp>
    typename LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::cend() const
    {
        return m_data.cbegin() + m_size;
    }

    // Element access --------------------------------------------------------------

    template <class T, int N, class Cmp>
    const T& LimitedPriorityQueue<T, N, Cmp>::top() const
    {
        return m_data[0];
    }

    template <class T, int N, class Cmp>
    const T& LimitedPriorityQueue<T, N, Cmp>::bottom() const
    {
        return m_data[m_size - 1];
    }

    template <class T, int N, class Cmp>
    T& LimitedPriorityQueue<T, N, Cmp>::top()
    {
        return m_data[0];
    }

    template <class T, int N, class Cmp>
    T& LimitedPriorityQueue<T, N, Cmp>::bottom()
    {
        return m_data[m_size - 1];
    }

    // Capacity --------------------------------------------------------------------

    template <class T, int N, class Cmp>
    bool LimitedPriorityQueue<T, N, Cmp>::empty() const
    {
        return m_size == 0;
    }

    template <class T, int N, class Cmp>
    bool LimitedPriorityQueue<T, N, Cmp>::full() const
    {
        return m_size == capacity();
    }

    template <class T, int N, class Cmp>
    size_t LimitedPriorityQueue<T, N, Cmp>::size() const
    {
        return m_size;
    }

    template <class T, int N, class Cmp>
    size_t LimitedPriorityQueue<T, N, Cmp>::capacity() const
    {
        return m_capacity;
    }

    // Modifiers -------------------------------------------------------------------
    template <class T, int N, class Cmp>
    bool LimitedPriorityQueue<T, N, Cmp>::pushImpl(const T& _value, T** _addr)
    {
        if (empty())
        {
            if (capacity() > 0)
            {
                *_addr = &m_data.front();
                ++m_size;
                return true;
            }
            return false; // No capacity to insert
        }
        iterator it = internal::upperBound(begin(), end(), _value, m_comp);
        if (it == end())
        {
            if (!full())
            {
                *_addr = &(*it);
                ++m_size;
                return true;
            }
            return false;
        }
        if (full())
        {
            internal::copyBackward(it, end() - 1, end());
            *_addr = &(*it);
        }
        else
        {
            internal::copyBackward(it, end(), end() + 1);
            *_addr = &(*it);
            ++m_size;
        }
        return true;
    }

    template <class T, int N, class Cmp>
    bool LimitedPriorityQueue<T, N, Cmp>::push(T&& _value)
    {
        T* addr;
        if (pushImpl(_value, &addr)) // Needs to insert
        {
            *addr = std::forward<T>(_value);
            return true;
        }
        return false;
    }

    template <class T, int N, class Cmp>
    bool LimitedPriorityQueue<T, N, Cmp>::push(const T& _value)
    {
        T* addr;
        if (pushImpl(_value, &addr)) // Needs to insert
        {
            *addr = _value;
            return true;
        }
        return false;
    }

    template <class T, int N, class Cmp>
    void LimitedPriorityQueue<T, N, Cmp>::pop()
    {
        --m_size;
    }

    template <class T, int N, class Cmp>
    void LimitedPriorityQueue<T, N, Cmp>::reserve(const int _capacity)
    {
        PONCA_ASSERT(_capacity >= 0);
        PONCA_ASSERT(_capacity <= N);
        m_capacity = _capacity;
        if (m_size > _capacity)
            m_size = _capacity;
    }

    template <class T, int N, class Cmp>
    void LimitedPriorityQueue<T, N, Cmp>::clear()
    {
        m_size = 0;
    }

    // Data ------------------------------------------------------------------------

    template <class T, int N, class Cmp>
    const typename LimitedPriorityQueue<T, N, Cmp>::container_type& LimitedPriorityQueue<T, N, Cmp>::container() const
    {
        return m_data;
    }
} // namespace Ponca
