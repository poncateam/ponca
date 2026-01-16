/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeKNearestIterator.h"

namespace Ponca {

/*!
 * \brief Extension of the Query class that allows to read the result of a k-nearest neighbor search on the KdTree.
 *
 *  Ouput result of a `KdTreeBase::kNearestNeighbors` query request.
 *
 *  \see KdTreeBase
 */
template <typename Traits,
          template <typename, typename> typename IteratorType,
          typename QueryType>
class KdTreeKNearestQueryBase : public KdTreeQuery<Traits>, public QueryType
{
public:
    using DataPoint      = typename Traits::DataPoint;
    using IndexType      = typename Traits::IndexType;
    using Scalar         = typename DataPoint::Scalar;
    using VectorType     = typename DataPoint::VectorType;
    using QueryAccelType = KdTreeQuery<Traits>;
    using Iterator       = IteratorType<typename Traits::IndexType, typename Traits::DataPoint>;
    using Self           = KdTreeKNearestQueryBase<Traits, IteratorType, QueryType>;

    inline KdTreeKNearestQueryBase(const KdTreeBase<Traits>* kdtree, IndexType k, typename QueryType::InputType input) :
            KdTreeQuery<Traits>(kdtree), QueryType(k, input) { }

    /// \brief Call the k-nearest neighbors query with new input and neighbor number parameters.
    inline Self& operator()(typename QueryType::InputType input, IndexType k)
    {
        return QueryType::template operator()<Self>(input, k);
    }
    /// \brief Call the k-nearest neighbors query with new input parameter.
    inline Self& operator()(typename QueryType::InputType input)
    {
        return QueryType::template operator()<Self>(input);
    }

    /// \brief Returns an iterator to the beginning of the k-nearest neighbors query.
    inline Iterator begin(){
        QueryAccelType::reset();
        QueryType::reset();
        this->search();
        return Iterator(QueryType::m_queue.begin());
    }

    /// \brief Returns an iterator to the end of the k-nearest neighbors query.
    inline Iterator end(){
        return Iterator(QueryType::m_queue.end());
    }

protected:
    inline void search(){
        KdTreeQuery<Traits>::search_internal(QueryType::getInputPosition(QueryAccelType::m_kdtree->points()),
                                             [](IndexType, IndexType){},
                                             [this](){return QueryType::descentDistanceThreshold();},
                                             [this](IndexType idx){return QueryType::skipIndexFunctor(idx);},
                                             [this](IndexType idx, IndexType, Scalar d){QueryType::m_queue.push({idx, d}); return false;}
        );
    }
};
/*!
 * \copybrief KdTreeKNearestQueryBase
 *
 * Ouput result of a `KdTreeBase::kNearestNeighbors` query made with the **index** of the point to evaluate.
 * \see KNearestIndexQuery
 */
template <typename Traits>
using KdTreeKNearestIndexQuery = KdTreeKNearestQueryBase< Traits, KdTreeKNearestIterator,
                                 KNearestIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>>;
/*!
 * \copybrief KdTreeKNearestQueryBase
 *
 * Ouput result of a `KdTreeBase::kNearestNeighbors` query made with the **position** of the point to evaluate.
 * \see KNearestPointQuery
 */
template <typename Traits>
using KdTreeKNearestPointQuery = KdTreeKNearestQueryBase< Traits, KdTreeKNearestIterator,
                                 KNearestPointQuery<typename Traits::IndexType, typename Traits::DataPoint>>;
} // namespace ponca
