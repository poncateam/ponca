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

namespace Ponca {

//!
//! \brief The limited_priority_queue class is similar to std::priority_queue
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
template<class T,
         class CompareT = std::less<T>>
class limited_priority_queue
{
public:
    using value_type      = T;
    using container_type  = std::vector<T>;
    using compare         = CompareT;
    using iterator        = typename container_type::iterator;
    using const_iterator  = typename container_type::const_iterator;
    using this_type       = limited_priority_queue<T,CompareT>;

    // limited_priority_queue --------------------------------------------------
public:
    PONCA_MULTIARCH inline limited_priority_queue();
    PONCA_MULTIARCH inline limited_priority_queue(const this_type& other);
    PONCA_MULTIARCH inline explicit limited_priority_queue(int capacity);
    template<class InputIt>
    PONCA_MULTIARCH inline limited_priority_queue(int capacity, InputIt first, InputIt last);

    PONCA_MULTIARCH inline ~limited_priority_queue();

    PONCA_MULTIARCH inline limited_priority_queue& operator=(const this_type& other);

    // Iterator ----------------------------------------------------------------
public:
    PONCA_MULTIARCH inline iterator       begin();
    PONCA_MULTIARCH inline const_iterator begin() const;
    PONCA_MULTIARCH inline const_iterator cbegin() const;

    PONCA_MULTIARCH inline iterator       end();
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
    PONCA_MULTIARCH inline bool empty() const;
    PONCA_MULTIARCH inline bool full() const;
    PONCA_MULTIARCH inline size_t  size() const;
    PONCA_MULTIARCH inline size_t  capacity() const;

    // Modifiers ---------------------------------------------------------------
public:
    PONCA_MULTIARCH inline bool push(const T& value);
    PONCA_MULTIARCH inline bool push(T&& value);

    PONCA_MULTIARCH inline void pop();

    PONCA_MULTIARCH inline void reserve(int capacity);

    PONCA_MULTIARCH inline void clear();

    // Data --------------------------------------------------------------------
public:
    PONCA_MULTIARCH inline const container_type& container() const;

protected:
    container_type  m_c;
    compare         m_comp;
    size_t          m_size {0};
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// limited_priority_queue ------------------------------------------------------

template<class T, class Cmp>
limited_priority_queue<T,Cmp>::limited_priority_queue() :
    m_c(),
    m_comp(),
    m_size(0)
{
}


template<class T, class Cmp>
limited_priority_queue<T,Cmp>::limited_priority_queue(const this_type& other) :
    m_c(other.m_c),
    m_comp(other.m_comp),
    m_size(other.m_size)
{
}

template<class T, class Cmp>
limited_priority_queue<T,Cmp>::limited_priority_queue(int capacity) :
    m_c(capacity),
    m_comp(),
    m_size(0)
{
}

template<class T, class Cmp>
template<class InputIt>
limited_priority_queue<T,Cmp>::limited_priority_queue(int capacity, InputIt first, InputIt last) :
    m_c(capacity),
    m_comp(),
    m_size(0)
{
    for(InputIt it=first; it<last; ++it)
    {
        push(*it);
    }
}

template<class T, class Cmp>
limited_priority_queue<T,Cmp>::~limited_priority_queue()
{
}

template<class T, class Cmp>
limited_priority_queue<T,Cmp>& limited_priority_queue<T,Cmp>::operator=(const this_type& other)
{
    m_c    = other.m_c;
    m_comp = other.m_comp;
    m_size = other.m_size;
    return *this;
}

// Iterator --------------------------------------------------------------------

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::iterator limited_priority_queue<T,Cmp>::begin()
{
    return m_c.begin();
}

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::const_iterator limited_priority_queue<T,Cmp>::begin() const
{
    return m_c.begin();
}

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::const_iterator limited_priority_queue<T,Cmp>::cbegin() const
{
    return m_c.cbegin();
}

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::iterator limited_priority_queue<T,Cmp>::end()
{
    return m_c.begin() + m_size;
}

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::const_iterator limited_priority_queue<T,Cmp>::end() const
{
    return m_c.begin() + m_size;
}

template<class T, class Cmp>
typename limited_priority_queue<T,Cmp>::const_iterator limited_priority_queue<T,Cmp>::cend() const
{
    return m_c.cbegin() + m_size;
}

// Element access --------------------------------------------------------------

template<class T, class Cmp>
const T& limited_priority_queue<T,Cmp>::top() const
{
    return m_c[0];
}

template<class T, class Cmp>
const T& limited_priority_queue<T,Cmp>::bottom() const
{
    return m_c[m_size-1];
}

template<class T, class Cmp>
T& limited_priority_queue<T,Cmp>::top()
{
    return m_c[0];
}

template<class T, class Cmp>
T& limited_priority_queue<T,Cmp>::bottom()
{
    return m_c[m_size-1];
}

// Capacity --------------------------------------------------------------------

template<class T, class Cmp>
bool limited_priority_queue<T,Cmp>::empty() const
{
    return m_size == 0;
}

template<class T, class Cmp>
bool limited_priority_queue<T,Cmp>::full() const
{
    return m_size == capacity();
}

template<class T, class Cmp>
size_t limited_priority_queue<T,Cmp>::size() const
{
    return m_size;
}

template<class T, class Cmp>
size_t limited_priority_queue<T,Cmp>::capacity() const
{
    return m_c.size();
}

// Modifiers -------------------------------------------------------------------

template<class T, class Cmp>
bool limited_priority_queue<T,Cmp>::push(const T& value)
{
    if(empty())
    {
        if(capacity()>0)
        {
            m_c.front() = value;
            ++m_size;
            return true;
        }
    }
    else
    {
        iterator it = std::upper_bound(begin(), end(), value, m_comp);
        if(it==end())
        {
            if(!full())
            {
                *it = value;
                ++m_size;
                return true;
            }
        }
        else
        {
            if(full())
            {
                std::copy_backward(it, end()-1, end());
                *it = value;
            }
            else
            {
                std::copy_backward(it, end(), end()+1);
                *it = value;
                ++m_size;
            }
            return true;
        }
    }
    return false;
}

template<class T, class Cmp>
bool limited_priority_queue<T,Cmp>::push(T&& value)
{
    if(empty())
    {
        if(capacity()>0)
        {
            m_c.front() = std::move(value);
            ++m_size;
            return true;
        }
    }
    else
    {
        iterator it = std::upper_bound(begin(), end(), std::move(value), m_comp);
        if(it==end())
        {
            if(!full())
            {
                *it = std::move(value);
                ++m_size;
                return true;
            }
        }
        else
        {
            if(full())
            {
                std::copy_backward(it, end()-1, end());
                *it = std::move(value);
            }
            else
            {
                std::copy_backward(it, end(), end()+1);
                *it = std::move(value);
                ++m_size;
            }
            return true;
        }
    }
    return false;
}

template<class T, class Cmp>
void limited_priority_queue<T,Cmp>::pop()
{
    --m_size;
}

template<class T, class Cmp>
void limited_priority_queue<T,Cmp>::reserve(int capacity)
{
    if(m_size>capacity)
    {
        m_size = capacity;
    }
    m_c.resize(capacity);
}

template<class T, class Cmp>
void limited_priority_queue<T,Cmp>::clear()
{
    m_size = 0;
}

// Data ------------------------------------------------------------------------

template<class T, class Cmp>
const typename limited_priority_queue<T,Cmp>::container_type& limited_priority_queue<T,Cmp>::container() const
{
    return m_c;
}

}   
