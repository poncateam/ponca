/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template<typename Index>
class KdTreeNearestIterator
{
public:
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using value_type = Index;
    using pointer    = Index*;
    using reference  = Index&;

    inline KdTreeNearestIterator() = default;
    inline KdTreeNearestIterator(Index index) : m_index(index) {}
    virtual inline ~KdTreeNearestIterator() = default;

public:
    inline bool operator ==(const KdTreeNearestIterator& other) const
    {return m_index == other.m_index;}
    inline bool operator !=(const KdTreeNearestIterator& other) const
    {return m_index != other.m_index;}
    /// Postfix increment
    inline KdTreeNearestIterator operator++(int) {
        KdTreeNearestIterator tmp = *this;
        ++m_index;
        return tmp;
    }
    /// Prefix increment
    inline KdTreeNearestIterator& operator ++() {++m_index; return *this;}
    inline reference operator *() const {
        return const_cast<reference>(m_index);
    }

protected:
    Index m_index {-1};
};

} // namespace ponca
