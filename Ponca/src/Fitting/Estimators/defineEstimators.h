#pragma once

#include "poncaKernels.h"

#include "defineEnum.h"
#include "estimator.h"

namespace Estimators {


template <typename WeightFunc>
using Fit_ASO = Ponca::BasketDiff<
            Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
            Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using Estimator_ASO = Estimator<Fit_ASO<WeightFunc>, true>;

}

/*

// template <typename WeightFunc>
// using fit_PCA = Ponca::BasketDiff<
//             Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>,
//             Ponca::DiffType::FitSpaceDer,
//             Ponca::CovariancePlaneDer,
//             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using fit_PCA = Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>;

template <typename WeightFunc>
using fit_MeanPLANE = Ponca::BasketDiff<
            Ponca::Basket<PPAdapter, WeightFunc, Ponca::MeanPlaneFit>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::MeanPlaneDer,
            Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

// template <typename WeightFunc>
// using fit_MeanPLANE = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MeanPlaneFit>;

/// Spheres with ASO approaches
// template <typename WeightFunc>
// using fit_APSS = Ponca::BasketDiff<
//             Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>,
//             Ponca::DiffType::FitSpaceDer,
//             Ponca::OrientedSphereDer,
//             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using fit_APSS = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleOrientedSphereFit>;

// template <typename WeightFunc>
// using fit_Sphere = Ponca::BasketDiff<
//             Ponca::Basket<PPAdapter, WeightFunc, Ponca::SphereFit>,
//             Ponca::DiffType::FitSpaceDer,
//             Ponca::SphereFitDer,
//             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using fit_Sphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleSphereFit>;

// template <typename WeightFunc>
// using fit_UnorientedSphere = Ponca::BasketDiff<
//             Ponca::Basket<PPAdapter, WeightFunc, Ponca::UnorientedSphereFit>,
//             Ponca::DiffType::FitSpaceDer,
//             Ponca::UnorientedSphereDer,
//             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using fit_UnorientedSphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleUnorientedSphereFit>;

template <typename WeightFunc>
using fit_MeanCurvate = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MeanCurvatureFit>;

template <typename WeightFunc>
using fit_Cov2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::Covariance2DFit>;

template <typename WeightFunc>
using fit_NormCov2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::NormalCovariance2DFit>;

template <typename WeightFunc>
using fit_NormCov3D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::NormalCovariance3DFit>;

template<typename WeightFunc>
using fit_ShapeOperator = Ponca::Basket<PPAdapter, WeightFunc, Ponca::ShapeOperator2DFit>;


template <typename WeightFunc>
using fit_Ellipsoid = Ponca::BasketDiff<
            Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedEllipsoidFit>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

template <typename WeightFunc>
using fit_BOCylinder = Ponca::Basket<PPAdapter, WeightFunc, Ponca::BaseOrientedParabolicCylinderFit>;

template <typename WeightFunc>
using fit_FOCylinder = Ponca::Basket<PPAdapter, WeightFunc, Ponca::FullyOrientedParabolicCylinderFit>;

template <typename WeightFunc>
using fit_BCylinder = Ponca::Basket<PPAdapter, WeightFunc, Ponca::BaseParabolicCylinderFit>;

template <typename WeightFunc>
using fit_B2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::BaseEllipsoid2DFit>;

template <typename WeightFunc>
using fit_BO2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::BaseOrientedEllipsoid2DFit>;

template <typename WeightFunc>
using fit_FO2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::FullyOrientedEllipsoid2DFit>;

template <typename WeightFunc>
using fit_WaveJets = Ponca::Basket<PPAdapter, WeightFunc, Ponca::WaveJetsFit>;

template <typename WeightFunc>
using fit_OrientedWaveJets = Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedWaveJetsFit>;

using fit_CNC = Ponca::Basket<PPAdapter, ConstWeightFunc, Ponca::TriangleGeneration>;

using fit_Varifolds = Ponca::Basket<PPAdapter, VarifoldWeightFunc, Ponca::VarifoldsCovPlane>;

using fit_VarifoldsMeanPlane = Ponca::Basket<PPAdapter, VarifoldWeightFunc, Ponca::VarifoldsMeanPlane>;

template <typename WeightFunc>
using fit_Monge = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MongePatchFit>;

template <typename WeightFunc>
using fit_OrientedMonge = Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedMongePatchFit>;

template <typename WeightFunc>
using fit_quadric = Ponca::Basket<PPAdapter, WeightFunc, Ponca::QuadricFit>;
*/
