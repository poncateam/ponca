#pragma once

#include <Ponca/Fitting>
#include <Ponca/src/Fitting/gls.h>
#include "Ponca/src/Common/pointTypes.h"

/*!
  \file fittingPCH.h
  \brief Define the commonly used Fitting types for 3D.
*/

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
