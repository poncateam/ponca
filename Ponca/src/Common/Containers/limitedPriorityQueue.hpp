/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble
*/

namespace Ponca
{
    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>& LimitedPriorityQueue<T, N, Cmp>::operator=(const Self& other) =
        default;

    // Iterator --------------------------------------------------------------------

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::iterator LimitedPriorityQueue<T, N, Cmp>::begin()
    {
        return m_data.begin();
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::begin() const
    {
        return m_data.begin();
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::cbegin() const
    {
        return m_data.cbegin();
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::iterator LimitedPriorityQueue<T, N, Cmp>::end()
    {
        return m_data.begin() + m_size;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH typename LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::end()
        const
    {
        return m_data.begin() + m_size;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::const_iterator LimitedPriorityQueue<T, N, Cmp>::cend() const
    {
        return m_data.cbegin() + m_size;
    }

    // Element access --------------------------------------------------------------

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH const T& LimitedPriorityQueue<T, N, Cmp>::top() const
    {
        return m_data[0];
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH const T& LimitedPriorityQueue<T, N, Cmp>::bottom() const
    {
        return m_data[m_size - 1];
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH T& LimitedPriorityQueue<T, N, Cmp>::top()
    {
        return m_data[0];
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH T& LimitedPriorityQueue<T, N, Cmp>::bottom()
    {
        return m_data[m_size - 1];
    }

    // Capacity --------------------------------------------------------------------

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH [[nodiscard]] bool LimitedPriorityQueue<T, N, Cmp>::empty() const
    {
        return m_size == 0;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH [[nodiscard]] bool LimitedPriorityQueue<T, N, Cmp>::full() const
    {
        return m_size == capacity();
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH [[nodiscard]] size_t LimitedPriorityQueue<T, N, Cmp>::size() const
    {
        return m_size;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH [[nodiscard]] size_t LimitedPriorityQueue<T, N, Cmp>::capacity() const
    {
        return m_capacity;
    }

    // Modifiers -------------------------------------------------------------------
    template <class T, int N, class Cmp>
    PONCA_MULTIARCH bool LimitedPriorityQueue<T, N, Cmp>::push(T&& value)
    {
        if (empty())
        {
            if (capacity() > 0)
            {
                m_data.front() = std::forward<T>(value);
                ++m_size;
                return true;
            }
        }
        else
        {
            iterator it = Ponca::internal::upperBound(begin(), end(), value, m_comp);
            if (it == end())
            {
                if (!full())
                {
                    *it = std::forward<T>(value);
                    ++m_size;
                    return true;
                }
            }
            else
            {
                if (full())
                {
                    Ponca::internal::copyBackward(it, end() - 1, end());
                    *it = std::forward<T>(value);
                }
                else
                {
                    Ponca::internal::copyBackward(it, end(), end() + 1);
                    *it = std::forward<T>(value);
                    ++m_size;
                }
                return true;
            }
        }
        return false;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH void LimitedPriorityQueue<T, N, Cmp>::pop()
    {
        --m_size;
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH void LimitedPriorityQueue<T, N, Cmp>::reserve(const int capacity)
    {
        PONCA_ASSERT(capacity >= 0);
        PONCA_ASSERT(capacity <= N);
        m_capacity = capacity;
        if (m_size > capacity)
        {
            m_size = capacity;
        }
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH void LimitedPriorityQueue<T, N, Cmp>::clear()
    {
        m_size = 0;
    }

    // Data ------------------------------------------------------------------------

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH const LimitedPriorityQueue<T, N, Cmp>::container_type& LimitedPriorityQueue<T, N, Cmp>::container()
        const
    {
        return m_data;
    }
}
