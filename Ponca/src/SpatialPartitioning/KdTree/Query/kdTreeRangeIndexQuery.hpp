/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<typename Traits>
auto KdTreeRangeIndexQuery<Traits>::begin() -> Iterator
{
    QueryAccelType::reset();
    QueryType::reset();
    Iterator it(this);
    this->advance(it);
    return it;
}

template<typename Traits>
auto KdTreeRangeIndexQuery<Traits>::end() -> Iterator
{
    return Iterator(this, QueryAccelType::m_kdtree->point_count());
}

template<typename Traits>
void KdTreeRangeIndexQuery<Traits>::advance(Iterator& it)
{
    const auto& points  = QueryAccelType::m_kdtree->point_data();
    const auto& indices = QueryAccelType::m_kdtree->index_data();
    const auto& point   = points[QueryType::input()].pos();

    if (points.empty() || indices.empty())
        throw std::invalid_argument("Empty KdTree");

    for(IndexType i=it.m_start; i<it.m_end; ++i)
    {
        IndexType idx = indices[i];
        if(idx == QueryType::input()) continue;

        Scalar d = (point - points[idx].pos()).squaredNorm();
        if(d < QueryType::descentDistanceThreshold())
        {
            it.m_index = idx;
            it.m_start = i+1;
            return;
        }
    }

    if (KdTreeQuery<Traits>::search_internal(point,
                                             [&it](IndexType start, IndexType end)
                                             {
                                                 it.m_start = start;
                                                 it.m_end   = end;
                                             },
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [this](IndexType idx){return QueryType::input() == idx;},
                                             [&it](IndexType idx, IndexType i, Scalar)
                                             {
                                                 it.m_index = idx;
                                                 it.m_start = i+1;
                                                 return true;
                                             }
        ))
        it.m_index = static_cast<IndexType>(points.size());
}
