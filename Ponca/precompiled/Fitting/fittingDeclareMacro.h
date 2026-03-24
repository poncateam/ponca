/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file fittingDeclareMacro.h
 * \brief Make macros to declare the commonly used Fitting types for 3D point clouds
 */

#pragma once

#include <Ponca/Fitting>
#include <Ponca/src/Common/pointTypes.h>

template <typename Scalar>
using Point = Ponca::PointPositionNormal<Scalar, 3>;

template <typename Scalar>
using WeightSmoothFuncL = Ponca::DistWeightFunc<Point<Scalar>, Ponca::SmoothWeightKernel<Scalar>>;
template <typename Scalar>
using WeightConstantFuncL = Ponca::DistWeightFunc<Point<Scalar>, Ponca::ConstantWeightKernel<Scalar>>;
template <typename Scalar>
using NoWeightFuncG = Ponca::NoWeightFuncGlobal<Point<Scalar>>;
template <typename Scalar>
using NoWeightFuncL = Ponca::NoWeightFunc<Point<Scalar>>;

#define DEF_FIT(SCALAR, WEIGHT, ...) EXTERN template class Ponca::Basket<Point<SCALAR>, WEIGHT<SCALAR>, __VA_ARGS__>;

#define DEF_FIT_FOR_ALL_FUNCS(SCALAR, ...)            \
    DEF_FIT(SCALAR, WeightSmoothFuncL, __VA_ARGS__)   \
    DEF_FIT(SCALAR, WeightConstantFuncL, __VA_ARGS__) \
    DEF_FIT(SCALAR, NoWeightFuncL, __VA_ARGS__)       \
    DEF_FIT(SCALAR, NoWeightFuncG, __VA_ARGS__)

#define DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, ...) \
    DEF_FIT(SCALAR, WeightSmoothFuncL, __VA_ARGS__) \
    DEF_FIT(SCALAR, WeightConstantFuncL, __VA_ARGS__)

/* Fits that can't be precompiled because they lack some methods but are still used in testing. */
// #define INCOMPLETE_FITTING_DEF(SCALAR) \
//     DEF_FIT_FOR_ALL_FUNCS(SCALAR, Ponca::MeanPosition) \
//     DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::MeanPosition, Ponca::CovarianceFitBase) \
//     DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::CovarianceLineFit)

// #define FITTING_DIFF_DEF(SCALAR) \
//     EXTERN template class Ponca::BasketDiff< \
//         Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::GLSParam>, \
//         Ponca::FitScaleSpaceDer, Ponca::OrientedSphereDer, Ponca::GLSDer \
//     >; \
//     EXTERN template class Ponca::BasketDiff< \
//         Ponca::Basket<Point, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>, \
//         Ponca::FitScaleSpaceDer, Ponca::UnorientedSphereDer, Ponca::CurvatureEstimatorDer,  \
//         Ponca::NormalDerivativeWeingartenEstimator, Ponca::WeingartenCurvatureEstimatorDer \
//     >;

/* Common fit types to precompile */
#define FITTING_DEF(SCALAR)                                                           \
    DEF_FIT_FOR_ALL_FUNCS(SCALAR, Ponca::MeanPlaneFit)                                \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::Plane)                              \
    /* CovariancePlane fits */                                                        \
    DEF_FIT_FOR_ALL_FUNCS(SCALAR, Ponca::CovariancePlaneFit)                          \
    /* MongePatch fits */                                                             \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::MongePatchQuadraticFit)             \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::MongePatchRestrictedQuadraticFit)   \
    /* OrientedSphere fits */                                                         \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::OrientedSphereFit)                  \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::OrientedSphereFit, Ponca::GLSParam) \
    /* UnorientedSphere fits */                                                       \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::UnorientedSphereFit)                \
    DEF_FIT_FOR_ALL_WEIGHTED_FUNCS(SCALAR, Ponca::UnorientedSphereFit, Ponca::GLSParam)
