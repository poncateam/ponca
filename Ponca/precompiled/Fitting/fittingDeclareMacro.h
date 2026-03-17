#pragma once

#include <Ponca/Fitting>
#include <Ponca/src/Fitting/gls.h>
#include "Ponca/src/Common/pointTypes.h"

/*!
  \file fittingDeclareMacro.h
  \brief Make macros to declare the commonly used Fitting types for 3D point clouds
*/

#ifndef EXTERN
#    define EXTERN
#endif

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


/* Fits that can't be precompiled because they lack some methods but are still used in testing. */
#define INCOMPLETE_FITTING_DEF(SCALAR) \
    typedef Ponca::Basket<Point, WeightConstantFuncLocal, Ponca::MeanPosition> FitConstantLocal; \
    typedef Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>, Ponca::MeanPosition> FitNoWeightLocal; \
    typedef Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>, Ponca::MeanPosition> FitNoWeightGlobal;

/* Common fit types to precompile */
#define FITTING_DEF(SCALAR) \
    /* Mean */                                                                                            \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MeanPlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MeanPlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::MeanPlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::MeanPlaneFit>; \
    /* Covariance-based fits */                                                                                 \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::CovariancePlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::CovariancePlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::CovariancePlaneFit>; \
    EXTERN template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::CovariancePlaneFit>; \
    EXTERN template class Ponca::BasketDiff<                                                 \
        Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::CovariancePlaneFit>,  \
        Ponca::FitSpaceDer, Ponca::CovariancePlaneDer,                                       \
        Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator             \
    >;                                                                                       \
    /* MongePatch fits */                                                                                                       \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchQuadraticFit>;             \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchQuadraticFit>;             \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchRestrictedQuadraticFit>;   \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchRestrictedQuadraticFit>;   \
    /* OrientedSphereFit */                                                                                                     \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit>;                    \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit>;                  \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::GLSParam>;   \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::GLSParam>; \
    /* UnorientedSphereFit */                                                                                                     \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit>;                    \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::UnorientedSphereFit>;                  \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>;   \
    EXTERN template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>; \
    EXTERN template class Ponca::BasketDiff<                                                                               \
        Ponca::Basket< Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam >,            \
        Ponca::FitSpaceDer, Ponca::OrientedSphereDer, Ponca::GLSDer,                                                       \
        Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator, Ponca::WeingartenCurvatureEstimatorDer   \
    >;
