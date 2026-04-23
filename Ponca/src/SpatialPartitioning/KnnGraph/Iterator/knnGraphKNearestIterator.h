/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../indexSquaredDistance.h"
#include <cstddef>

namespace Ponca
{

    /*!
     *  \brief Input iterator to read the `KnnGraphKNearestQueryBase` object.
     *
     *  As this is an input iterator, we don't guarantee anything other than reading the values with it.
     *  If you need to operate on the values of this iterator with algorithms that relies on forward iterator
     * functionalities, you should copy the index values in an STL-like container.
     *
     *  \warning This iterator object should never be duplicated, as it is a proxy that holds a reference to the actual
     * data : The copy of this iterator would still point to the same KnnGraph reference. So, if the increment operator
     * is used on the iterator, the duplicate will also have its state updated. If we then call the increment operator
     * on the duplicate, the result will be an incorrect value.
     *
     *  \see KnnGraphKNearestQueryBase
     */
    template <typename Container, typename Index>
    class KnnGraphKNearestIterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Index;
        using pointer           = Index*;
        using reference         = const Index&;

        PONCA_MULTIARCH KnnGraphKNearestIterator(const Container* data, Index i) : m_data(data), m_i(i) {}

        /// \brief Inequality operand
        PONCA_MULTIARCH bool operator!=(const KnnGraphKNearestIterator& other) const { return m_i != other.m_i; }

        /// \brief Equality operand
        PONCA_MULTIARCH bool operator==(const KnnGraphKNearestIterator& other) const { return m_i == other.m_i; }

        /// \brief Equality operand
        PONCA_MULTIARCH KnnGraphKNearestIterator& operator++()
        {
            ++m_i;
            return *this;
        }

        /// \brief Dereference operator
        PONCA_MULTIARCH reference operator*() const { return (*m_data)[m_i]; }

        /// \brief Postfix increment
        PONCA_MULTIARCH inline KnnGraphKNearestIterator operator++(Index)
        {
            KnnGraphKNearestIterator tmp = *this;
            ++m_i;
            return tmp;
        }

        /// \brief Value increment
        PONCA_MULTIARCH inline void operator+=(const Index i) { m_i += i; }

    protected:
        const Container* m_data{nullptr};
        Index m_i{0};
    };
} // namespace Ponca
