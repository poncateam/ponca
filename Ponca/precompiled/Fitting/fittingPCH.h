#pragma once

#define PONCA_PCH

#include <Ponca/Fitting>
#include "tests/common/testUtils.h"

using namespace Ponca;

typedef PointPositionNormal<float, 3> Point;

typedef DistWeightFunc<Point, SmoothWeightKernel<float>> WeightSmoothFunc;
typedef DistWeightFunc<Point, ConstantWeightKernel<float>> WeightConstantFuncLocal;
typedef NoWeightFuncGlobal<Point> NoWeightFuncGlobalP;
typedef NoWeightFunc<Point> NoWeightFuncP;

typedef Basket<Point, WeightConstantFuncLocal, MeanPosition> FitConstantLocal;
typedef Basket<Point, NoWeightFuncP, MeanPosition> FitNoWeightLocal;
typedef Basket<Point, NoWeightFuncGlobalP, MeanPosition> FitNoWeightGlobal;

typedef Basket<Point, WeightSmoothFunc, CovariancePlaneFit> CovFitSmooth;
typedef Basket<Point, WeightConstantFuncLocal, CovariancePlaneFit> CovFitConstant;
typedef Basket<Point, NoWeightFuncP, CovariancePlaneFit> CovFitConstantNoWeight;
typedef Basket<Point, NoWeightFuncGlobalP, CovariancePlaneFit> CovFitConstantGlobal;

typedef Basket<Point, WeightSmoothFunc, MeanPlaneFit> MeanFitSmooth;
typedef Basket<Point, WeightConstantFuncLocal, MeanPlaneFit> MeanFitConstant;
typedef Basket<Point, NoWeightFuncP, MeanPlaneFit> MeanFitConstant2;
typedef Basket<Point, NoWeightFuncGlobalP, MeanPlaneFit> MeanFitConstantGlobal;

typedef Basket<Point, WeightSmoothFunc, CovariancePlaneFit, MongePatch> CovFitSmoothMongePatch;
typedef Basket<Point, WeightConstantFuncLocal, CovariancePlaneFit, MongePatch> CovFitConstantMongePatch;

//! [PlaneFitDerTypes]
using TestPlane = Basket<Point, WeightSmoothFunc, CovariancePlaneFit>;
using PlaneScaleDiff = BasketDiff<TestPlane, FitScaleDer, CovariancePlaneDer>;
using PlaneSpaceDiff = BasketDiff<TestPlane, FitSpaceDer, CovariancePlaneDer>;
using PlaneScaleSpaceDiff = BasketDiff<TestPlane, FitScaleSpaceDer, CovariancePlaneDer>;
//! [PlaneFitDerTypes]
