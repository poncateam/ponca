/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file Ponca/src/Common/Containers/iteratorUtils.h
 * \brief Provides functions on iterator that works (Removes CUDA call Host function in device warning)
 */
#pragma once

#include "../defines.h"

namespace Ponca::internal
{
    /*! \brief Searches for the first element in the partitioned range [first, last) which is ordered after value.
     *
     * \tparam ForwardIt Must meet the requirements of LegacyForwardIterator
     * \tparam Compare Must meet the requirements of BinaryPredicate. It is not required to satisfy Compare.
     * \param first Iterator defining the beginning of the range of elements to examine
     * \param last Iterator defining the end of the range of elements to examine
     * \param value Value to compare the elements to
     * \param comp binary predicate which returns True if the first argument is ordered before the second.
     * \return Iterator to the first element of the range [first, last) ordered after value, or last if no such element
     * is found.
     *
     * \see https://en.cppreference.com/cpp/algorithm/copy_backward
     */
    template <class ForwardIt, class T = typename std::iterator_traits<ForwardIt>::value_type, class Compare>
    PONCA_MULTIARCH ForwardIt upperBound(ForwardIt first, ForwardIt last, const T& value, Compare comp)
    {
        ForwardIt it;
        typename std::iterator_traits<ForwardIt>::difference_type count, step;
        count = last - first;

        while (count > 0)
        {
            it   = first;
            step = count / 2;
            it += step;

            if (!comp(value, *it))
            {
                first = ++it;
                count -= step + 1;
            }
            else
                count = step;
        }

        return first;
    }

    /*! \brief Copies the elements from the range [first, last) to another range ending at d_last. The elements are
     * copied in reverse order (the last element is copied first), but their relative order is preserved.
     *
     * \warning The behavior is undefined if d_last is within (first, last].
     *
     * \tparam BidirIt1 must meet the requirements of LegacyBidirectionalIterator
     * \tparam BidirIt2 must meet the requirements of LegacyBidirectionalIterator
     * \param first Iterator defining the beginning of the range of elements to copy from (the source)
     * \param first Iterator defining the end of the range of elements to copy from (the source)
     * \param d_last The end of the destination range
     *
     * \see https://en.cppreference.com/cpp/algorithm/copy_backward
     */
    template <class BidirIt1, class BidirIt2>
    PONCA_MULTIARCH BidirIt2 copyBackward(BidirIt1 first, BidirIt1 last, BidirIt2 d_last)
    {
        while (first != last)
            *(--d_last) = *(--last);
        return d_last;
    }
} // namespace Ponca::internal