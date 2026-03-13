#include "fittingPCH.h"

#define FIT(SCALAR) \
    /* Mean */                                                                                     \
    extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MeanPosition>; \
    extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MeanPosition>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::MeanPosition>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::MeanPosition>; \
                                                                                                   \
    extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MeanPlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MeanPlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::MeanPlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::MeanPlaneFit>; \
    /* Covariance-based fits */                                                                          \
    extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::CovariancePlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::CovariancePlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::CovariancePlaneFit>; \
    extern template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::CovariancePlaneFit>; \
    // extern template class BasketDiff<                                                                                           \
    //     Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::CovariancePlaneFit>,                                     \
    //     Ponca::FitSpaceDer, Ponca::CovariancePlaneDer, Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator \
    // >;                                                                                                                          \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchQuadraticFit>;             \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchQuadraticFit>;             \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchRestrictedQuadraticFit>;   \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchRestrictedQuadraticFit>;   \
    // /* OrientedSphereFit */                                                                                                                 \
    // extern template class Basket<Point<SCALAR>, WeightSmoothFuncl<SCALAR>, Ponca::OrientedSphereFit, Ponca::OrientedSphereSpaceDer>;        \
    // extern template class Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::OrientedSphereSpaceDer>;      \
    // extern template class Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::OrientedSphereScaleSpaceDer>;   \
    // extern template class Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::OrientedSphereScaleSpaceDer>; \
    //                                                                                                                                         \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::OrientedSphereFit, Ponca::OrientedSphereSpaceDer, Ponca::MlsSphereFitDer>;      \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, OrientedSphereSpaceDer, Ponca::MlsSphereFitDer>;             \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::OrientedSphereFit, Ponca::OrientedSphereScaleSpaceDer, Ponca::MlsSphereFitDer>; \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::OrientedSphereScaleSpaceDer, Ponca::MlsSphereFitDer>; \
    // /* UnorientedSphereFit*/                                                                                                       \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>;    \
    // extern template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>;  \
    // extern template class Ponca::BasketDiff<                                                                                       \
    //     Ponca::Basket<MyPoint<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>,                    \
    //     Ponca::FitSpaceDer, Ponca::OrientedSphereDer, Ponca::GLSDer,                                                               \
    //     Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator, Ponca::WeingartenCurvatureEstimatorDer           \
    // >;

FIT(float)
FIT(double)
FIT(long double)

#undef FIT
