/**
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Thibault Lejemble - Main source code
 \author Nicolas Mellado   - Improved error handling and bound check
*/

#pragma once

#include <array>

#include "../defines.h" //STD_SAFE_AT

namespace Ponca
{

    /// Stack with fixed-size storage
    ///
    /// \warning This class does not destroy elements that are removed from the stack, but rather when
    /// they are overriden by another element (e.g. when running pop() and then push(). This makes it
    /// not appropriate to store smart pointers.
    ///
    template <class T, int N>
    class Stack
    {
    public:
        /// Type of value stored in the Stack
        using ValueType = T;

        PONCA_MULTIARCH inline Stack();

        /// Read access to the top element of the Stack
        PONCA_MULTIARCH inline const T& top() const;
        /// Write access to the top element of the Stack
        PONCA_MULTIARCH inline T& top();

        /// Is the stack empty
        PONCA_MULTIARCH inline bool empty() const;
        /// Get the number of elements in the Stack
        PONCA_MULTIARCH inline int size() const;

        /// Add an element on top of the stack.
        /// \throw std::out_of_range when the Stack is full, only if compiled with PONCA_DEBUG
        PONCA_MULTIARCH inline void push(const T& value);

        /// Add an element with default initialization
        PONCA_MULTIARCH inline void push();
        /// Pop the last element of the Stack
        /// \note In practice the element is not destroyed, only the size of the Stack is reduced. This
        /// makes the Stack not appropriate to store `shared_ptr`: counters will not be decremented.
        PONCA_MULTIARCH inline void pop();
        /// Clear the stack content.
        /// \note Shares the same limitation than Stack::pop()
        PONCA_MULTIARCH inline void clear();

    protected:
        /// Number of elements in the Stack
        int m_size;
        /// Fixed-size data buffer
        std::array<T, N> m_data;
    };
} // namespace Ponca

#include "stack.hpp"
