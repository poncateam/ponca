/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble - Main source code
 \author Nicolas Mellado   - Improved error handling and bound check
*/

#pragma once

#include <array>

namespace Ponca {

/// Stack with fixed-size storage
///
/// \warning This class does not destroy elements that are removed from the stack, but rather when
/// they are overriden by another element (e.g. when running pop() and then push(). This makes it
/// not appropriate to store smart pointers.
///
///
/// \ingroup common
template<class T, int N>
class Stack
{
public:
    /// Type of value stored in the Stack
    using ValueType = T;

    inline Stack();

    /// Read access to the top element of the Stack
    inline const T& top() const;
    /// Write access to the top element of the Stack
    inline       T& top();

    /// Is the stack empty
    inline bool empty() const;
    /// Get the number of elements in the Stack
    inline int  size() const;

    /// Add an element on top of the stack.
    /// \throw std::out_of_range when the Stack is full, only if not compiled with NDEBUG
    inline void push(const T& value);

    /// Add an element with default initialization
    inline void push();
    /// Pop the last element of the Stack
    /// \note In practice the element is not destroyed, only the size of the Stack is reduced. This
    /// makes the Stack not appropriate to store `shared_ptr`: counters will not be decremented.
    inline void pop();
    /// Clear the stack content.
    /// \note Shares the same limitation than Stack::pop()
    inline void clear();

protected:
    /// Number of elements in the Stack
    int m_size;
    /// Fixed-size data buffer
    std::array<T,N> m_data;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class T, int N>
Stack<T,N>::Stack() :
    m_size(0),
    m_data()
{
}

template<class T, int N>
const T& Stack<T,N>::top() const
{
    return STD_SAFE_AT(m_data,m_size-1);
}

template<class T, int N>
T& Stack<T,N>::top()
{
    return STD_SAFE_AT(m_data,m_size-1);
}

template<class T, int N>
bool Stack<T,N>::empty() const
{
    return m_size==0;
}

template<class T, int N>
int Stack<T,N>::size() const
{
    return m_size;
}

template<class T, int N>
void Stack<T,N>::push(const T& value)
{
    STD_SAFE_AT(m_data,m_size) = value;
    ++m_size;
}

template<class T, int N>
void Stack<T,N>::push()
{
    ++m_size;
}

template<class T, int N>
void Stack<T,N>::pop()
{
    --m_size;
}

template<class T, int N>
void Stack<T,N>::clear()
{
    m_size = 0;
}

} // namespace Ponca

