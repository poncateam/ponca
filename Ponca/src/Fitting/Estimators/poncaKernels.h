#pragma once

#include "baseType.h"

#include <Eigen/Dense>
#include <Ponca/Common>
#include <Ponca/Fitting>


namespace Ponca::Estimators {
    using SmoothWeightFunc      = Ponca::DistWeightFunc<PPAdapter, Ponca::SmoothWeightKernel<Scalar> >;
    // using SixSmoothWeightFunc   = Ponca::DistWeightFunc<PPAdapter, Ponca::SixSmoothWeightKernel<Scalar> >;
    // using FourSmoothWeightFunc  = Ponca::DistWeightFunc<PPAdapter, Ponca::FourSmoothWeightKernel<Scalar> >;
    // using ThreeSmoothWeightFunc = Ponca::DistWeightFunc<PPAdapter, Ponca::ThreeSmoothWeightKernel<Scalar> >;
    // using ConstWeightFunc       = Ponca::DistWeightFunc<PPAdapter, Ponca::ConstantWeightKernel<Scalar> >;
    // using WendlandWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::WendlandWeightKernel<Scalar> >;
    // using SingularWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::SingularWeightKernel<Scalar> >;
    // using VarifoldWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::VarifoldWeightKernel<Scalar> >;
    // using ExponentialWeightFunc = Ponca::DistWeightFunc<PPAdapter, Ponca::ExponentialWeightKernel<Scalar> >;
    // using RationalWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::RationalWeightKernel<Scalar> >;
}