#pragma once

#include <limits>
#include "./defines.h"

namespace Ponca {

/// \brief Associates an index with a distance
/// \ingroup spatialpartitioning
struct IndexSquaredDistance
{
    //// Index of the closest point
    int index {-1};

    /// Distance to the closest point
    SPScalar squared_distance { std::numeric_limits<SPScalar>::max() };

    /// Comparison operator based on squared_distance
    inline bool operator < (const IndexSquaredDistance& other) const
    { return squared_distance < other.squared_distance; }
};

}   
