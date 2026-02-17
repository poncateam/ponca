/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include PONCA_MULTIARCH_INCLUDE_STD(iterator)

namespace Ponca{
    /*!
      \brief ComputeObject is a virtual object that represents an algorithm which can be used with the compute functions.

      The compute(begin, end) and computeWithIds(ids, points) methods must be implemented by the inheriting class.
      \note The compute(container) that is defined in this structure can be reused in the inheriting class by adding
      "using ComputeObject<Self>::compute;" to make it accessible
 */
    template <typename Derived>
    struct ComputeObject {
    protected:
        /// \brief Retrieve the top layer object
        /// Returns a reference to the derived class so that we can use its overwritten methods
        PONCA_MULTIARCH Derived& derived() { return static_cast<Derived&>(*this); }
    public:

#ifdef PONCA_CPU_ARCH
        /*! \brief Convenience function for STL-like container
         *
         * \note This method is only accessible when using a CPU architecture (PONCA_CPU_ARCH = true)
         * \tparam Container And STL-Like container
         * \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
         */
        template <typename Container>
        FIT_RESULT compute(const Container& c) {
            return derived().compute(std::begin(c), std::end(c));
        }
#endif

        /*! \brief Convenience function for STL-like iterators
            \tparam IteratorBegin The beginning of the iterator (std::begin(iterator)
            \tparam IteratorEnd   The end of the iterator (std::end(iterator)
        */
        template <typename IteratorBegin, typename IteratorEnd>
        PONCA_MULTIARCH inline FIT_RESULT compute(const IteratorBegin& /*begin*/, const IteratorEnd& /*end*/) {
            return UNDEFINED;
        };

        /*! \brief Convenience function to iterate over a subset of samples in a PointContainer
            \tparam IndexRange STL-Like range storing indices of the neighbors
            \tparam PointContainer STL-like container storing the points
            \see #compute(const IteratorBegin& begin, const IteratorEnd& end)
        */
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIds(IndexRange /*ids*/, const PointContainer& /*points*/) {
            return UNDEFINED;
        };
    }; // struct ComputeObject

}

