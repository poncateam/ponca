/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <typename Traits>
auto KdTreeNearestPointQuery<Traits>::begin()
    -> KdTreeNearestIterator<IndexType>
{
    QueryAccelType::reset();
    QueryType::reset();
    this->search();
    return KdTreeNearestIterator<IndexType>(QueryType::m_nearest);
}

template <typename Traits>
auto KdTreeNearestPointQuery<Traits>::end()
    -> KdTreeNearestIterator<IndexType>
{
    return KdTreeNearestIterator<IndexType>(QueryType::m_nearest + 1);
}

template <typename Traits>
void KdTreeNearestPointQuery<Traits>::search()
{
    KdTreeQuery<Traits>::search_internal(QueryType::input(),
                                         [](IndexType, IndexType){},
                                         [this](){return QueryType::descentDistanceThreshold();},
                                         [](IndexType){return false;},
                                         [this](IndexType idx, IndexType, Scalar d)
                                         {
                                             QueryType::m_nearest = idx;
                                             QueryType::m_squared_distance = d;
                                             return false;
                                         }
    );
}
