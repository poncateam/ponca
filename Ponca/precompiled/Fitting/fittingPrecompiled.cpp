/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file Ponca/precompiled/Fitting/fittingPrecompiled.cpp
 * \brief Precompiled target that defines the commonly used Fitting types for 3D point clouds.
 */

#include "fittingPCH.h"

#define EXTERN
#include "fittingDeclareMacro.h"

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef FITTING_DEF
#undef EXTERN

#define CALL_FIT(SCALAR, WEIGHT, ...)                                  \
    {                                                                  \
        Ponca::Basket<Point<SCALAR>, WEIGHT<SCALAR>, __VA_ARGS__> fit; \
        fit.setNeighborFilter({point});                                \
        fit.compute(points);                                           \
    }

#define CALL_FIT_FOR_ALL_FUNCS(SCALAR, ...)            \
    CALL_FIT(SCALAR, WeightSmoothFuncL, __VA_ARGS__)   \
    CALL_FIT(SCALAR, WeightConstantFuncL, __VA_ARGS__) \
    CALL_FIT(SCALAR, NoWeightFuncL, __VA_ARGS__)       \
    CALL_FIT(SCALAR, NoWeightFuncG, __VA_ARGS__)

#define CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, ...) \
    CALL_FIT(SCALAR, WeightSmoothFuncL, __VA_ARGS__) \
    CALL_FIT(SCALAR, WeightConstantFuncL, __VA_ARGS__)

/* Common fit types to precompile */
#define FITTING_CALL(SCALAR)                                                           \
    CALL_FIT_FOR_ALL_FUNCS(SCALAR, Ponca::MeanPlaneFit)                                \
    /* CovariancePlane fits */                                                         \
    CALL_FIT_FOR_ALL_FUNCS(SCALAR, Ponca::CovariancePlaneFit)                          \
    /* MongePatch fits */                                                              \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::MongePatchQuadraticFit)             \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::MongePatchRestrictedQuadraticFit)   \
    /* OrientedSphere fits */                                                          \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::OrientedSphereFit)                  \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::OrientedSphereFit, Ponca::GLSParam) \
    /* UnorientedSphere fits */                                                        \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::UnorientedSphereFit)                \
    CALL_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::UnorientedSphereFit, Ponca::GLSParam)

namespace Ponca::internal
{
    void main(const int /*argc*/, char** /*argv*/)
    {
        // Explicitly call all the functions to pre-instantiate them
        {
            const Point<float> point;
            std::vector<Point<float>> points = {point, point};
            FITTING_CALL(float)
        }
        {
            const Point<double> point;
            std::vector<Point<double>> points = {point, point};
            FITTING_CALL(double)
        }
        {
            const Point<long double> point;
            std::vector<Point<long double>> points = {point, point};
            FITTING_CALL(long double)
        }
    }
} // namespace Ponca::internal
