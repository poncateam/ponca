#include "fittingPCH.h"

/* Fits that can't be precompiled because they lack some methods but are still used in testing. */
#define INCOMPLETE_FITTING_DEF(SCALAR) \
    typedef Ponca::Basket<Point, WeightConstantFuncLocal, Ponca::MeanPosition> FitConstantLocal; \
    typedef Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>, Ponca::MeanPosition> FitNoWeightLocal; \
    typedef Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>, Ponca::MeanPosition> FitNoWeightGlobal;

/* Common fit types to precompile */
#define FITTING_DEF(SCALAR) \
    /* Mean */                                                                                     \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MeanPlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MeanPlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::MeanPlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::MeanPlaneFit>; \
    /* Covariance-based fits */                                                                          \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::CovariancePlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::CovariancePlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, NoWeightFuncL<SCALAR>      , Ponca::CovariancePlaneFit>; \
    template class Ponca::Basket<Point<SCALAR>, NoWeightFuncG<SCALAR>      , Ponca::CovariancePlaneFit>; \
    template class Ponca::BasketDiff<                                                        \
        Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::CovariancePlaneFit>,  \
        Ponca::FitSpaceDer, Ponca::CovariancePlaneDer,                                       \
        Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator             \
    >;                                                                                       \
    /* MongePatch fits */                                                                                                \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchQuadraticFit>;             \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchQuadraticFit>;             \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>  , Ponca::MongePatchRestrictedQuadraticFit>;   \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::MongePatchRestrictedQuadraticFit>;   \
    /* OrientedSphereFit */                                                                                              \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit>;                    \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit>;                  \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::GLSParam>;   \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::OrientedSphereFit, Ponca::GLSParam>; \
    /* UnorientedSphereFit */                                                                                              \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit>;                    \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::UnorientedSphereFit>;                  \
    template class Ponca::Basket<Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>;   \
    template class Ponca::Basket<Point<SCALAR>, WeightConstantFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam>; \
    template class Ponca::BasketDiff<                                                                                      \
        Ponca::Basket< Point<SCALAR>, WeightSmoothFuncL<SCALAR>, Ponca::UnorientedSphereFit, Ponca::GLSParam >,            \
        Ponca::FitSpaceDer, Ponca::OrientedSphereDer, Ponca::GLSDer,                                                       \
        Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator, Ponca::WeingartenCurvatureEstimatorDer   \
    >;

FITTING_DEF(float)
FITTING_DEF(double)
FITTING_DEF(long double)

#undef FITTING_DEF
