/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "enums.h"
#include <iostream>

namespace Ponca
{
    template < class DataPoint, class _WFunctor, typename T>
    class MLS : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES

        protected:
        int mlsIter = 1;
        VectorType lastPosMLS = VectorType::Zero();
    public:
        PONCA_EXPLICIT_CAST_OPERATORS(MLS,mls)

        /// Sets the number of mls iterations that will be done during the compute process.
        inline void setIterMLS(const int _maxIter) {
            mlsIter = _maxIter;
        }

        /// Returns the last projected position of the current eval position that was processed in the compute method
        inline VectorType getLastMLSPosition() const {
            return lastPosMLS;
        }

        /// Needs ti have a weight func set before calling this method
        /// Returns the projected position after the mls iteration
        template <typename IndexRange, typename PointContainer>
        PONCA_MULTIARCH inline FIT_RESULT computeWithIdsMLS(IndexRange ids, const PointContainer& points);
    };
#include "mls.hpp"
}