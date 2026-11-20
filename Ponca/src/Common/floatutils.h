/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "defines.h"
#include <Eigen/Core>
#include <concepts>

#ifndef PONCA_ABS_EPS
#    define PONCA_ABS_EPS Scalar(1)
#endif

#ifndef PONCA_REL_EPS
#    define PONCA_REL_EPS Eigen::NumTraits<Scalar>::dummy_precision()
#endif

namespace Ponca
{
    namespace internal
    {
        template <typename Derived, typename OtherDerived, typename Scalar = typename Derived::RealScalar>
        PONCA_MULTIARCH inline bool isApprox(const Eigen::DenseBase<Derived>& x,
                                             const Eigen::DenseBase<OtherDerived>& y, Scalar abs = PONCA_ABS_EPS,
                                             Scalar rel = PONCA_REL_EPS)
        {
            // From eigen isApprox_selector (Core/Fuzzy.h:25 for current Eigen's version)
            typename Eigen::internal::nested_eval<Derived, 3>::type X(x.derived());
            typename Eigen::internal::nested_eval<OtherDerived, 3>::type Y(y.derived());

            const Scalar eps2 = rel * rel;
            const Scalar diff = (X - Y).cwiseAbs2().sum();

            // Add an absolute test here to avoid instabilities near 0
            if (diff <= rel * abs)
                return true;

            return diff <= eps2 * X.cwiseAbs2().sum() && diff <= eps2 * Y.cwiseAbs2().sum();
        }

        template <std::floating_point Scalar>
        PONCA_MULTIARCH inline bool isApprox(Scalar x, Scalar y, Scalar abs = PONCA_ABS_EPS, Scalar rel = PONCA_REL_EPS)
        {
            const Scalar eps2 = rel * rel;
            const Scalar diff = (x - y) * (x - y);

            if (diff <= eps2 * abs * abs)
                return true;

            return diff <= eps2 * x * x && diff <= eps2 * y * y;
        }

        template <std::integral Scalar>
        PONCA_MULTIARCH inline bool isApprox(Scalar x, Scalar y, Scalar abs = Scalar(0), Scalar rel = Scalar(0))
        {
            return x == y;
        }
    } // namespace internal
} // namespace Ponca
