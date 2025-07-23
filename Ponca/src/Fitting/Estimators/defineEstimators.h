#pragma once

#include "poncaKernels.h"
#include "estimator.h"
#include "baseType.h"

namespace Ponca::Estimators {
    /// Distance to PCA plane : Mean curvature
    template <typename WeightFunc>
    using Fit_PCA = Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>;

    /// Distance to PCA plane : Curvature Tensor0
    // template <typename WeightFunc>
    // using Fit_PCA = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::CovariancePlaneDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Point Set Surfaces : Curvature Tensor
    template <typename WeightFunc>
    using Fit_PSS = Ponca::BasketDiff<
            Fit_PCA<WeightFunc>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::CovariancePlaneDer,
    Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Algebraic Point Set Surfaces : Mean curvature
    template <typename WeightFunc>
    using Fit_APSS = Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>;

    /// Algebraic Point Set Surfaces : Curvature Tensor
    // template <typename WeightFunc>
    // using Fit_APSS = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::OrientedSphereDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;


    /// Algebraic Shape Operator : Curvature Tensor
    template <typename WeightFunc>
    using Fit_ASO = Ponca::BasketDiff<
                Fit_APSS<WeightFunc>,
                Ponca::DiffType::FitSpaceDer,
                Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
                Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

    /// Sphere
    // template <typename WeightFunc>
    // using Fit_Sphere = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::SphereFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::SphereFitDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    template <typename WeightFunc>
    using Fit_Sphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SphereFit>;

}
