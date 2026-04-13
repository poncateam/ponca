/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#ifdef _PONCA_COMPILE_DEFINITION
#   define _PRECOMPILED_PONCA_EXTERN
#else
#   define _PRECOMPILED_PONCA_EXTERN extern
#endif

#include "../Ponca"


namespace Ponca
{
    // Examples for ponca-basic cpu
    _PRECOMPILED_PONCA_EXTERN template class Basket<PointPositionNormal<double, 3>, DistWeightFunc<PointPositionNormal<double, 3>, SmoothWeightKernel<double>>, OrientedSphereFit, GLSParam>;;
    _PRECOMPILED_PONCA_EXTERN template class Basket<PointPositionNormal<double, 3>, DistWeightFunc<PointPositionNormal<double, 3>, SmoothWeightKernel<double>>, UnorientedSphereFit, GLSParam>;
    _PRECOMPILED_PONCA_EXTERN template class BasketDiff<Basket<PointPositionNormal<double, 3>, DistWeightFunc<PointPositionNormal<double, 3>, SmoothWeightKernel<double>>, OrientedSphereFit, GLSParam>, FitSpaceDer, OrientedSphereDer, GLSDer, CurvatureEstimatorDer, NormalDerivativeWeingartenEstimator, WeingartenCurvatureEstimatorDer>;
    _PRECOMPILED_PONCA_EXTERN template class Basket<PointPositionNormal<double, 3>, DistWeightFunc<PointPositionNormal<double, 3>, SmoothWeightKernel<double>>, SphereFit, GLSParam>;
};