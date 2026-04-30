/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../defines.h"
#include "../enums.h"

namespace Ponca
{
    /**
     * \brief SingleEvaluationScheme
     *
     * Simple evaluation scheme that calls compute on the given fit
     */
    struct SingleEvaluationScheme
    {
    public:
#ifdef PONCA_CPU_ARCH
        /*! \brief Convenience function for STL-like container
         *
         * \note This method is only accessible when using a CPU architecture (PONCA_CPU_ARCH = true)
         * \tparam ComputeObject Compute object to call the function on
         * \tparam Container And STL-Like container
         * \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
         */
        template <typename ComputeObject, typename Container>
        FIT_RESULT compute(ComputeObject& co, const Container& c) const
        {
            return co.compute(std::begin(c), std::end(c));
        }
#endif

        /*! \brief Convenience function for STL-like iterators
            \tparam ComputeObject Compute object to call the function on
            \tparam IteratorBegin The beginning of the iterator (std::begin(iterator)
            \tparam IteratorEnd   The end of the iterator (std::end(iterator)
        */
        //! [SingleEvaluationScheme Compute Definition]
        template <typename ComputeObject, typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline FIT_RESULT compute(ComputeObject& co, const IteratorBegin& begin,
                                                  const IteratorEnd& end) const
        {
            return co.compute(begin, end);
        };
        //! [SingleEvaluationScheme Compute Definition]

        /*! \brief Convenience function to iterate over a subset of samples in a PointContainer
            \tparam ComputeObject Compute object to call the function on
            \tparam IndexRange STL-Like range storing indices of the neighbors
            \tparam PointContainer STL-like container storing the points
            \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
        */
        template <typename ComputeObject, typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(ComputeObject& co, IndexRange ids,
                                                         const PointContainer& points) const
        {
            return co.computeWithIds(ids, points);
        };

    private:
    };
} // namespace Ponca
