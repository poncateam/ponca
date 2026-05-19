/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../indexSquaredDistance.h"

namespace Ponca
{

    /*!
     *  \brief Base iterator class for NeighborGraphOneConnectedQuery.
     *
     *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
     *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator
     * functionalities, you should copy the index values in an STL-like container.
     *
     *  \warning This iterator object should never be duplicated, as it is a proxy that holds a reference to the actual
     * data : the copy of this iterator would still point to the same NeighborGraph reference. So, modifying one will
     * modify the others and can result in incorrect values.
     */
    template <typename ContainerPtr, typename Index>
    class NeighborGraphOneConnectedIterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Index;
        using pointer           = ContainerPtr;
        using reference         = const Index&;
        using Self              = NeighborGraphOneConnectedIterator<ContainerPtr, Index>;

        PONCA_MULTIARCH NeighborGraphOneConnectedIterator(ContainerPtr data) : m_data(data) {}

        PONCA_MULTIARCH NeighborGraphOneConnectedIterator(ContainerPtr data, Index i) : m_data(data), m_i(i) {}

        /// \brief Inequality operand
        PONCA_MULTIARCH bool operator!=(const NeighborGraphOneConnectedIterator& other) const
        {
            return m_i != other.m_i;
        }

        /// \brief Equality operand
        PONCA_MULTIARCH bool operator==(const NeighborGraphOneConnectedIterator& other) const
        {
            return m_i == other.m_i;
        }

        /// \brief Equality operand
        PONCA_MULTIARCH NeighborGraphOneConnectedIterator& operator++()
        {
            ++m_i;
            return *this;
        }

        /// \brief Dereference operator
        PONCA_MULTIARCH reference operator*() const { return m_data[m_i]; }

        /// \brief Postfix increment
        PONCA_MULTIARCH inline Self operator++(Index)
        {
            Self tmp = *this;
            ++m_i;
            return tmp;
        }

        /// \brief Value increment
        PONCA_MULTIARCH inline void operator+=(const Index i) { m_i += i; }

        /// \brief Plus operator
        PONCA_MULTIARCH inline Self operator+(const Index i) { return Self(m_data, m_i + i); }

    protected:
        ContainerPtr m_data;
        Index m_i{0};
    };
} // namespace Ponca
