/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "defines.h"
#include "defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(iterator)

namespace Ponca
{

template <class P, class W>
struct CNC : BasketBase {
    using random = Eigen::internal::random<int>;

    int nbTry = 1;

    template <typename IndexRange, typename PointContainer>
    FIT_RESULT computeWithIds(IndexRange ids, const PointContainer& points){
        float curvature = 0;
        int lengthIds = ids.size();

        for (int i = 0; i < nbTry; i++) {
            int t1 = points[ids[random(0, lengthIds)]];
            int t2 = points[ids[random(0, lengthIds)]];
            int t3 = points[ids[random(0, lengthIds)]];
            curvature += estimate(t1, t2, t3);
        }

        curvature /= nbTry;
    }

    template <typename PointContainer>
    PONCA_MULTIARCH inline
    FIT_RESULT compute(const PointContainer& points) {
        float curvature = 0;
        int lengthPoints = points.size();

        for (int i = 0; i < nbTry; i++) {
            int t1 = points[random(0, lengthPoints)];
            int t2 = points[random(0, lengthPoints)];
            int t2 = points[random(0, lengthPoints)];
            curvature += estimate(t1, t2, t3);
        }

        curvature /= nbTry;
    };
};
