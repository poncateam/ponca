/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble
*/

#pragma once

#include <cstddef>
#include <cassert>
#include <array>
#include <algorithm>
#include <functional>
#include <Ponca/src/Common/Containers/iteratorUtils.h>

#include "../defines.h"

namespace Ponca
{

    //!
    //! \brief The LimitedPriorityQueue class is similar to std::priority_queue
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
    //! \tparam T The data type stored in the queue
    //! \tparam N The maximum capacity of the queue
    //! \tparam CompareT A binary predicate used to sort the queue. Default to less
    template <class T, int N, class CompareT = std::less<T>>
    class LimitedPriorityQueue
    {
        static_assert(N > 0, "The capacity must be strictly positive");

    public:
        using value_type     = T;
        using container_type = std::array<T, N>;
        // using container_type = std::vector<T>;
        using compare        = CompareT;
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        using Base           = LimitedPriorityQueue<T, N, CompareT>;

        // LimitedPriorityQueue --------------------------------------------------
    public:
        PONCA_MULTIARCH inline LimitedPriorityQueue();
        PONCA_MULTIARCH inline LimitedPriorityQueue(const Base& other);
        PONCA_MULTIARCH inline explicit LimitedPriorityQueue(int capacity);
        template <class InputIt>
        PONCA_MULTIARCH inline LimitedPriorityQueue(int capacity, InputIt first, InputIt last);

        PONCA_MULTIARCH inline ~LimitedPriorityQueue();

        PONCA_MULTIARCH inline LimitedPriorityQueue& operator=(const Base& other);

        // Iterator ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH inline iterator begin();
        PONCA_MULTIARCH inline const_iterator begin() const;
        PONCA_MULTIARCH inline const_iterator cbegin() const;

        PONCA_MULTIARCH inline iterator end();
        PONCA_MULTIARCH inline const_iterator end() const;
        PONCA_MULTIARCH inline const_iterator cend() const;

        // Element access ----------------------------------------------------------
    public:
        PONCA_MULTIARCH inline const T& top() const;
        PONCA_MULTIARCH inline const T& bottom() const;

        PONCA_MULTIARCH inline T& top();
        PONCA_MULTIARCH inline T& bottom();

        // Capacity ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH [[nodiscard]] inline bool empty() const;
        PONCA_MULTIARCH [[nodiscard]] inline bool full() const;
        PONCA_MULTIARCH [[nodiscard]] inline size_t size() const;
        PONCA_MULTIARCH [[nodiscard]] inline size_t capacity() const;

        // Modifiers ---------------------------------------------------------------
    public:
        PONCA_MULTIARCH inline bool push(T&& value);
        PONCA_MULTIARCH inline void pop();
        PONCA_MULTIARCH inline void reserve(int capacity);
        PONCA_MULTIARCH inline void clear();

        // Data --------------------------------------------------------------------
    public:
        PONCA_MULTIARCH inline const container_type& container() const;

    protected:
        container_type m_data{};
        compare m_comp;
        size_t m_size{0};     //!< The current size of the Queue
        size_t m_capacity{0}; //!< The capacity of the limited queue
    };

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    // LimitedPriorityQueue ------------------------------------------------------

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue() : m_comp()
    {
        assert((m_capacity <= N));
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue(const Base& other)
        : m_data(other.m_data), m_comp(other.m_comp), m_size(other.m_size), m_capacity(other.m_capacity)
    {
        assert((m_capacity <= N));
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue(const int capacity)
        : m_comp(), m_capacity(capacity)
    {
        assert((capacity >= 0));
        assert((capacity <= N));
    }

    template <class T, int N, class Cmp>
    template <class InputIt>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue(const int capacity, InputIt first,
                                                                          InputIt last)
        : m_comp(), m_capacity(capacity)
    {
        for (InputIt it = first; it < last; ++it)
        {
            push(*it);
        }
        assert((capacity >= 0));
        assert((capacity <= N));
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::~LimitedPriorityQueue() = default;

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>& LimitedPriorityQueue<T, N, Cmp>::operator=(const Base& other) =
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
        assert(capacity >= 0);
        assert(capacity <= N);
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

} // namespace Ponca
