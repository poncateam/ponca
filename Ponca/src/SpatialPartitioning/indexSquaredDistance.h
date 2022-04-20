/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <limits>
#include "./defines.h"

namespace Ponca {

/// \brief Associates an index with a distance
/// \ingroup spatialpartitioning
template<typename Scalar>
struct IndexSquaredDistance
{
    //// Index of the closest point
    int index {-1};

    /// Distance to the closest point
    Scalar squared_distance { std::numeric_limits<Scalar>::max() };

    /// Comparison operator based on squared_distance
    inline bool operator < (const IndexSquaredDistance& other) const
    { return squared_distance < other.squared_distance; }
};

}   
