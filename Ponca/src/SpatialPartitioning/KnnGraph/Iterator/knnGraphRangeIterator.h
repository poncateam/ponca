/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

template <typename Traits>
class KnnGraphRangeQuery;

template <typename Traits>
class KnnGraphRangeIterator
{
protected:
    friend class KnnGraphRangeQuery<Traits>;

public:
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    using value_type = int;
    using pointer    = int*;
    using reference  = int&;

    inline KnnGraphRangeIterator(KnnGraphRangeQuery<Traits>* query, int index = -1) : m_query(query), m_index(index) {}

public:
    bool operator != (const KnnGraphRangeIterator& other) const{
        return m_index != other.m_index;
    }

    void operator ++ (){
        m_query->advance(*this);
    }

    inline reference operator *() const {
        return const_cast<reference>(m_index);
    }

protected:
    KnnGraphRangeQuery<Traits>* m_query {nullptr};
    int m_index {-1};
};

} // namespace Ponca
