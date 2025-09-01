#pragma once

#include "poncaKernels.h"
#include "estimator.h"

namespace Ponca::Estimators {
    /// Distance to PCA plane : Mean curvature
    template <typename WeightFunc, typename DataPoint>
    using Fit_PCA = Basket<DataPoint, WeightFunc, Ponca::CovariancePlaneFit>;

    /// Distance to PCA plane : Curvature Tensor
    // template <typename WeightFunc>
    // using Fit_PCA_Diff = Ponca::BasketDiff<
    //             Fit_PCA<WeightFunc>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::CovariancePlaneDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Point Set Surfaces : Curvature Tensor
    template <typename WeightFunc, typename DataPoint>
    using Fit_PSS = Ponca::BasketDiff<
            Fit_PCA<WeightFunc, DataPoint>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::CovariancePlaneDer,
    Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Algebraic Point Set Surfaces : Mean curvature
    template <typename WeightFunc, typename DataPoint>
    using Fit_APSS = Ponca::Basket<DataPoint, WeightFunc, Ponca::OrientedSphereFit>;

    /// Algebraic Point Set Surfaces : Curvature Tensor
    // template <typename WeightFunc>
    // using Fit_APSS_Diff = Ponca::BasketDiff<
    //             Fit_APSS<WeightFunc>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::OrientedSphereDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;


    /// Algebraic Shape Operator : Curvature Tensor
    template <typename WeightFunc, typename DataPoint>
    using Fit_ASO = Ponca::BasketDiff<
                Fit_APSS<WeightFunc, DataPoint>,
                Ponca::DiffType::FitSpaceDer,
                Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
                Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Sphere : Mean curvature
    template <typename WeightFunc, typename DataPoint>
    using Fit_Sphere = Ponca::Basket<DataPoint, WeightFunc, Ponca::SphereFit>;
    /// Sphere : Curvature Tensor
    // template <typename WeightFunc>
    // using Fit_Sphere_Diff = Ponca::BasketDiff<
    //             Fit_Sphere,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::SphereFitDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
}
