/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble
*/

#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include "../defines.h"

namespace Ponca
{

    //!
    //! \brief The limitedPriorityQueue class is similar to std::priority_queue
    //! but has a limited capacity and handles the comparison differently.
    //!
    //! In case the capacity is reached, the container is full and push() do not
    //! insert a new element if its priority is lower than the current minimal one.
    //!
    //! The comparison predicate must return true is the first argument has
    //! priority on the second one.
    //!
    //! The element with the highest priority is the last one that is pop out, but
    //! is the first one that is iterated through.
    //!
    //! Example 1:
    //! Using std::less as comparison predicate, we have the following situation:
    //!
    //!     full()      = false
    //!     empty()     = false
    //!     capacity()  = 6
    //!     size()      = 4
    //!     top()       = 1
    //!     bottom()    = 9
    //!
    //!     pop() removes the value 9
    //!     push(4) adds the value 4
    //!
    //!     begin            end
    //!       v               v
    //!     +---+---+---+---+---+---+
    //!     | 1 | 3 | 8 | 9 |   |   |
    //!     +---+---+---+---+---+---+
    //!       ^           ^
    //!      top        bottom
    //!
    //!
    //!
    //! Example 2:
    //! Using std::greater as comparison predicate, we have the following situation:
    //!
    //!     full()      = true
    //!     empty()     = false
    //!     capacity()  = 6
    //!     size()      = 6
    //!     top()       = 9
    //!     bottom()    = 2
    //!
    //!     begin                    end
    //!       v                       v
    //!     +---+---+---+---+---+---+
    //!     | 9 | 8 | 6 | 4 | 3 | 2 |
    //!     +---+---+---+---+---+---+
    //!       ^                   ^
    //!      top                bottom
    //!
    //!     pop() removes the value 2
    //!     push(5) adds the value 5 and remove the value 2
    //!     push(0) do nothing
    //!
    template <class T, class CompareT = std::less<T>>
    class limitedPriorityQueue
    {
    public:
        using value_type     = T;
        using container_type = std::vector<T>;
        using compare        = CompareT;
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        using this_type      = limitedPriorityQueue<T, CompareT>;

        // limitedPriorityQueue --------------------------------------------------
    public:
        PONCA_MULTIARCH_HOST inline limitedPriorityQueue();
        PONCA_MULTIARCH_HOST inline limitedPriorityQueue(const this_type& other);
        PONCA_MULTIARCH_HOST inline explicit limitedPriorityQueue(int capacity);
        template <class InputIt>
        PONCA_MULTIARCH_HOST inline limitedPriorityQueue(int capacity, InputIt first, InputIt last);

        PONCA_MULTIARCH_HOST inline ~limitedPriorityQueue();

        PONCA_MULTIARCH inline limitedPriorityQueue& operator=(const this_type& other);

        // Iterator ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH_HOST inline iterator begin();
        PONCA_MULTIARCH_HOST inline const_iterator begin() const;
        PONCA_MULTIARCH_HOST inline const_iterator cbegin() const;

        PONCA_MULTIARCH_HOST inline iterator end();
        PONCA_MULTIARCH_HOST inline const_iterator end() const;
        PONCA_MULTIARCH_HOST inline const_iterator cend() const;

        // Element access ----------------------------------------------------------
    public:
        PONCA_MULTIARCH_HOST inline const T& top() const;
        PONCA_MULTIARCH_HOST inline const T& bottom() const;

        PONCA_MULTIARCH_HOST inline T& top();
        PONCA_MULTIARCH_HOST inline T& bottom();

        // Capacity ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH [[nodiscard]] inline bool empty() const;
        PONCA_MULTIARCH_HOST [[nodiscard]] inline bool full() const;
        PONCA_MULTIARCH [[nodiscard]] inline size_t size() const;
        PONCA_MULTIARCH_HOST [[nodiscard]] inline size_t capacity() const;

        // Modifiers ---------------------------------------------------------------
    public:
        PONCA_MULTIARCH_HOST inline bool push(const T& value);
        PONCA_MULTIARCH_HOST inline bool push(T&& value);
        PONCA_MULTIARCH inline void pop();
        PONCA_MULTIARCH_HOST inline void reserve(int capacity);
        PONCA_MULTIARCH inline void clear();

        // Data --------------------------------------------------------------------
    public:
        PONCA_MULTIARCH inline const container_type& container() const;

    protected:
        container_type m_c;
        compare m_comp;
        size_t m_size{0};
    };

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    // limitedPriorityQueue ------------------------------------------------------

    template <class T, class Cmp>
    limitedPriorityQueue<T, Cmp>::limitedPriorityQueue() : m_c(), m_comp(), m_size(0)
    {
    }

    template <class T, class Cmp>
    limitedPriorityQueue<T, Cmp>::limitedPriorityQueue(const this_type& other)
        : m_c(other.m_c), m_comp(other.m_comp), m_size(other.m_size)
    {
    }

    template <class T, class Cmp>
    limitedPriorityQueue<T, Cmp>::limitedPriorityQueue(int capacity) : m_c(capacity), m_comp(), m_size(0)
    {
    }

    template <class T, class Cmp>
    template <class InputIt>
    limitedPriorityQueue<T, Cmp>::limitedPriorityQueue(int capacity, InputIt first, InputIt last)
        : m_c(capacity), m_comp(), m_size(0)
    {
        for (InputIt it = first; it < last; ++it)
        {
            push(*it);
        }
    }

    template <class T, class Cmp>
    limitedPriorityQueue<T, Cmp>::~limitedPriorityQueue()
    {
    }

    template <class T, class Cmp>
    limitedPriorityQueue<T, Cmp>& limitedPriorityQueue<T, Cmp>::operator=(const this_type& other)
    {
        m_c    = other.m_c;
        m_comp = other.m_comp;
        m_size = other.m_size;
        return *this;
    }

    // Iterator --------------------------------------------------------------------

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::iterator limitedPriorityQueue<T, Cmp>::begin()
    {
        return m_c.begin();
    }

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::const_iterator limitedPriorityQueue<T, Cmp>::begin() const
    {
        return m_c.begin();
    }

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::const_iterator limitedPriorityQueue<T, Cmp>::cbegin() const
    {
        return m_c.cbegin();
    }

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::iterator limitedPriorityQueue<T, Cmp>::end()
    {
        return m_c.begin() + m_size;
    }

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::const_iterator limitedPriorityQueue<T, Cmp>::end() const
    {
        return m_c.begin() + m_size;
    }

    template <class T, class Cmp>
    typename limitedPriorityQueue<T, Cmp>::const_iterator limitedPriorityQueue<T, Cmp>::cend() const
    {
        return m_c.cbegin() + m_size;
    }

    // Element access --------------------------------------------------------------

    template <class T, class Cmp>
    const T& limitedPriorityQueue<T, Cmp>::top() const
    {
        return m_c[0];
    }

    template <class T, class Cmp>
    const T& limitedPriorityQueue<T, Cmp>::bottom() const
    {
        return m_c[m_size - 1];
    }

    template <class T, class Cmp>
    T& limitedPriorityQueue<T, Cmp>::top()
    {
        return m_c[0];
    }

    template <class T, class Cmp>
    T& limitedPriorityQueue<T, Cmp>::bottom()
    {
        return m_c[m_size - 1];
    }

    // Capacity --------------------------------------------------------------------

    template <class T, class Cmp>
    bool limitedPriorityQueue<T, Cmp>::empty() const
    {
        return m_size == 0;
    }

    template <class T, class Cmp>
    bool limitedPriorityQueue<T, Cmp>::full() const
    {
        return m_size == capacity();
    }

    template <class T, class Cmp>
    size_t limitedPriorityQueue<T, Cmp>::size() const
    {
        return m_size;
    }

    template <class T, class Cmp>
    size_t limitedPriorityQueue<T, Cmp>::capacity() const
    {
        return m_c.size();
    }

    // Modifiers -------------------------------------------------------------------

    template <class T, class Cmp>
    bool limitedPriorityQueue<T, Cmp>::push(const T& value)
    {
        if (empty())
        {
            if (capacity() > 0)
            {
                m_c.front() = value;
                ++m_size;
                return true;
            }
        }
        else
        {
            iterator it = std::upper_bound(begin(), end(), value, m_comp);
            if (it == end())
            {
                if (!full())
                {
                    *it = value;
                    ++m_size;
                    return true;
                }
            }
            else
            {
                if (full())
                {
                    std::copy_backward(it, end() - 1, end());
                    *it = value;
                }
                else
                {
                    std::copy_backward(it, end(), end() + 1);
                    *it = value;
                    ++m_size;
                }
                return true;
            }
        }
        return false;
    }

    template <class T, class Cmp>
    bool limitedPriorityQueue<T, Cmp>::push(T&& value)
    {
        if (empty())
        {
            if (capacity() > 0)
            {
                m_c.front() = std::move(value);
                ++m_size;
                return true;
            }
        }
        else
        {
            iterator it = std::upper_bound(begin(), end(), std::move(value), m_comp);
            if (it == end())
            {
                if (!full())
                {
                    *it = std::move(value);
                    ++m_size;
                    return true;
                }
            }
            else
            {
                if (full())
                {
                    std::copy_backward(it, end() - 1, end());
                    *it = std::move(value);
                }
                else
                {
                    std::copy_backward(it, end(), end() + 1);
                    *it = std::move(value);
                    ++m_size;
                }
                return true;
            }
        }
        return false;
    }

    template <class T, class Cmp>
    void limitedPriorityQueue<T, Cmp>::pop()
    {
        --m_size;
    }

    template <class T, class Cmp>
    void limitedPriorityQueue<T, Cmp>::reserve(int capacity)
    {
        if (m_size > capacity)
        {
            m_size = capacity;
        }
        m_c.resize(capacity);
    }

    template <class T, class Cmp>
    void limitedPriorityQueue<T, Cmp>::clear()
    {
        m_size = 0;
    }

    // Data ------------------------------------------------------------------------

    template <class T, class Cmp>
    const typename limitedPriorityQueue<T, Cmp>::container_type& limitedPriorityQueue<T, Cmp>::container() const
    {
        return m_c;
    }

} // namespace Ponca
