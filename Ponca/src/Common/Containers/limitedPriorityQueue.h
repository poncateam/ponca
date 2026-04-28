/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble
*/

#pragma once

#include <cstddef>
#include <array>
#include <algorithm>
#include <functional>

#include "../defines.h"
#include "./iteratorUtils.h"
#include "../Assert.h"

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
        using compare        = CompareT;
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        using Self           = LimitedPriorityQueue<T, N, CompareT>;

        // LimitedPriorityQueue --------------------------------------------------
    public:
        PONCA_MULTIARCH LimitedPriorityQueue();
        PONCA_MULTIARCH LimitedPriorityQueue(const Self& other);
        PONCA_MULTIARCH explicit LimitedPriorityQueue(int capacity);
        template <class InputIt>
        PONCA_MULTIARCH LimitedPriorityQueue(int capacity, InputIt first, InputIt last);

        PONCA_MULTIARCH inline ~LimitedPriorityQueue() = default;

        PONCA_MULTIARCH LimitedPriorityQueue& operator=(const Self& other);

        // Iterator ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH iterator begin();
        PONCA_MULTIARCH const_iterator begin() const;
        PONCA_MULTIARCH const_iterator cbegin() const;

        PONCA_MULTIARCH iterator end();
        PONCA_MULTIARCH const_iterator end() const;
        PONCA_MULTIARCH const_iterator cend() const;

        // Element access ----------------------------------------------------------
    public:
        PONCA_MULTIARCH const T& top() const;
        PONCA_MULTIARCH const T& bottom() const;

        PONCA_MULTIARCH T& top();
        PONCA_MULTIARCH T& bottom();

        // Capacity ----------------------------------------------------------------
    public:
        PONCA_MULTIARCH [[nodiscard]] inline bool empty() const;
        PONCA_MULTIARCH [[nodiscard]] inline bool full() const;
        PONCA_MULTIARCH [[nodiscard]] inline size_t size() const;
        PONCA_MULTIARCH [[nodiscard]] inline size_t capacity() const;

        // Modifiers ---------------------------------------------------------------
    protected:
        PONCA_MULTIARCH inline bool pushImpl(const T& _value, T** _addr);

    public:
        PONCA_MULTIARCH bool push(T&& _value);
        PONCA_MULTIARCH bool push(const T& _value);
        PONCA_MULTIARCH void pop();
        PONCA_MULTIARCH void reserve(int _capacity);
        PONCA_MULTIARCH void clear();

        // Data --------------------------------------------------------------------
    public:
        PONCA_MULTIARCH const container_type& container() const;

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
        PONCA_ASSERT((m_capacity <= N));
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue(const Self& other)
        : m_data(other.m_data), m_comp(other.m_comp), m_size(other.m_size), m_capacity(other.m_capacity)
    {
        PONCA_ASSERT((m_capacity <= N));
    }

    template <class T, int N, class Cmp>
    PONCA_MULTIARCH LimitedPriorityQueue<T, N, Cmp>::LimitedPriorityQueue(const int capacity)
        : m_comp(), m_capacity(capacity)
    {
        PONCA_ASSERT((capacity >= 0));
        PONCA_ASSERT((capacity <= N));
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
        PONCA_ASSERT((capacity >= 0));
        PONCA_ASSERT((capacity <= N));
    }
} // namespace Ponca

#include "limitedPriorityQueue.hpp"
