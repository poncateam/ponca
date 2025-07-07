#pragma once

#include "poncaKernels.h"
#include "estimator.h"

namespace Estimators {
    /// Plane fit
    template <typename WeightFunc>
    using FitPlane = Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>;
    template <typename WeightFunc>
    using Fit_PSS = Ponca::BasketDiff<
            FitPlane<WeightFunc>,
            Ponca::DiffType::FitSpaceDer,
            Ponca::CovariancePlaneDer,
    Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    template <typename WeightFunc>
    using Estimator_PSS = Estimator<Fit_PSS<WeightFunc>, true>;


    /// APSS
    template <typename WeightFunc>
    using Fit_APSS = Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>;
    template <typename WeightFunc>
    using Estimator_APSS = Estimator<Fit_APSS<WeightFunc>, true>;
    /// Spheres with ASO approaches
    // template <typename WeightFunc>
    // using Fit_APSS = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::OrientedSphereDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;


    /// ASO fit
    template <typename WeightFunc>
    using Fit_ASO = Ponca::BasketDiff<
                Fit_APSS<WeightFunc>,
                Ponca::DiffType::FitSpaceDer,
                Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
                Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    template <typename WeightFunc>
    using Estimator_ASO = Estimator<Fit_ASO<WeightFunc>, true>;


    /// PCA
    // template <typename WeightFunc>
    // using Fit_PCA = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::CovariancePlaneDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    template <typename WeightFunc>
    using Fit_PCA = Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovariancePlaneFit>;
    template <typename WeightFunc>
    using Estimator_PCA = Estimator<Fit_PCA<WeightFunc>, true>;


    /// Sphere
    // template <typename WeightFunc>
    // using Fit_Sphere = Ponca::BasketDiff<
    //             Ponca::Basket<PPAdapter, WeightFunc, Ponca::SphereFit>,
    //             Ponca::DiffType::FitSpaceDer,
    //             Ponca::SphereFitDer,
    //             Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    template <typename WeightFunc>
    using Fit_Sphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SphereFit>;
    template <typename WeightFunc>
    using Estimator_Sphere = Estimator<Fit_Sphere<WeightFunc>, false>;

}
