#pragma once

#include "../Common/concepts.h"
#include "enums.h"

namespace Ponca
{
    template <typename T>
    concept HasLocalFrame = T::hasLocalFrame;

    template <typename T>
    concept Is3D = (T::Dim == 3);

    /// \brief This concept ensures that the default types and accessors in a Basket are well-formed
    ///
    /// This concept is implemented by BasketUnitBase, which is set by default in the Basket.
    template <typename T>
    concept ProvidesBasketUnitBase = ProvidesCommonTypes<T> && requires(T t, const T ct) {
        typename T::NeighborFilter;
        t.setNeighborFilter(typename T::NeighborFilter{});
        { ct.getNeighborFilter() } -> std::convertible_to<typename T::NeighborFilter>;

        { ct.getCurrentState() } -> std::same_as<FIT_RESULT>;

        { ct.getWeightSum() } -> std::same_as<typename T::Scalar>;
        { ct.isStable() } -> std::same_as<bool>;
        { ct.isReady() } -> std::same_as<bool>;
        { ct.getNumNeighbors() } -> std::integral;
    };

    /// \brief This concept ensures that the default types and accessors in a BasketDiff are well-formed
    ///
    /// This concept is implemented by BasketDiffUnitBase, which is set by default in the BasketDiff.
    template <typename T>
    concept ProvidesBasketDiffUnitBase = ProvidesBasketUnitBase<T> && requires(T t, const T ct) {
        { ct.isScaleDer() } -> std::same_as<bool>;
        { ct.isSpaceDer() } -> std::same_as<bool>;
        { ct.derDimension() } -> std::integral;
    };

    template <typename T>
    concept ProvidesProjectionOperator = requires(T t, const T ct, typename T::VectorType v, typename T::Scalar s) {
        ct.projectionOperator();

        { ct.projectionOperator().project(v) } -> std::same_as<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesImplicitPrimitive =
        ProvidesProjectionOperator<T> && requires(T t, const T ct, typename T::VectorType v, typename T::Scalar s) {
            t.implicitPrimitive();
            ct.implicitPrimitive();

            { ct.implicitPrimitive().potential() } -> std::same_as<typename T::Scalar>;
            { ct.implicitPrimitive().potential(v) } -> std::same_as<typename T::Scalar>;

            { ct.implicitPrimitive().primitiveGradient() } -> std::convertible_to<typename T::VectorType>;
            { ct.implicitPrimitive().primitiveGradient(v) } -> std::same_as<typename T::VectorType>;

            t.implicitPrimitive().changeBasis(v);
        };

    template <typename T>
    concept ProvidesAlgebraicSphere =
        ProvidesImplicitPrimitive<T> && requires(T t, const T ct, typename T::VectorType v, typename T::Scalar s) {
            ct.algebraicSphere();

            { ct.algebraicSphere().isPlane() } -> std::same_as<bool>;
            { ct.algebraicSphere().isValid() } -> std::same_as<bool>;
            { ct.algebraicSphere().isNormalized() } -> std::same_as<bool>;
            { ct.algebraicSphere().isApprox(ct, s) } -> std::same_as<bool>;

            { t.algebraicSphere().applyPrattNorm() } -> std::same_as<bool>;
            { ct.algebraicSphere().prattNorm() } -> std::same_as<typename T::Scalar>;
            { ct.algebraicSphere().prattNorm2() } -> std::same_as<typename T::Scalar>;

            { ct.algebraicSphere().radius() } -> std::convertible_to<typename T::Scalar>;
            { ct.algebraicSphere().center() } -> std::convertible_to<typename T::VectorType>;
        };

    template <typename T>
    concept ProvidesAlgebraicSphereDerivative = requires(const T ct) {
        ct.algebraicSphereDer();

        { ct.algebraicSphereDer().dPotential() } -> std::convertible_to<typename T::ScalarArray>;
    };

    template <typename T>
    concept ProvidesMeanPosition = requires(const T ct) {
        ct.meanPosition();

        { ct.meanPosition().barycenter() } -> std::convertible_to<typename T::VectorType>;
        { ct.meanPosition().barycenterDistance() } -> std::same_as<typename T::Scalar>;
        // This is only for subclasses
        // TODO: Make this function public / private ?
        // { ct.meanPosition().barycenterLocal() } -> std::convertible_to<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesMeanPositionDerivative = requires(const T ct) {
        ct.meanPositionDer();

        { ct.meanPositionDer().barycenterDerivatives() } -> std::convertible_to<typename T::VectorArray>;
    };

    template <typename T>
    concept ProvidesMeanNormal = requires(const T ct) {
        ct.meanNormal();

        { ct.meanNormal().meanNormalVector() } -> std::convertible_to<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesMeanNormalDer = requires(const T ct) {
        ct.meanNormalDer();

        { ct.meanNormal().dMeanNormal() } -> std::convertible_to<typename T::VectorArray>;
    };

    template <typename T>
    concept ProvidesGLSParam = requires(const T ct) {
        ct.glsParam();

        { ct.glsParam().tau() } -> std::same_as<typename T::Scalar>;
        { ct.glsParam().eta() } -> std::same_as<typename T::VectorType>;
        { ct.glsParam().kappa() } -> std::same_as<typename T::Scalar>;

        { ct.glsParam().tau_normalized() } -> std::same_as<typename T::Scalar>;
        { ct.glsParam().eta_normalized() } -> std::same_as<typename T::VectorType>;
        { ct.glsParam().kappa_normalized() } -> std::same_as<typename T::Scalar>;

        { ct.glsParam().fitness() } -> std::same_as<typename T::Scalar>;
        { ct.glsParam().compareTo(ct, true) } -> std::same_as<typename T::Scalar>;
    };

    template <typename T>
    concept ProvidesGLSDer = requires(const T ct) {
        ct.glsDer();

        { ct.glsDer().tau() } -> std::same_as<typename T::ScalarArray>;
        { ct.glsDer().eta() } -> std::same_as<typename T::VectorTypeArray>;
        { ct.glsDer().kappa() } -> std::same_as<typename T::ScalarArray>;

        { ct.glsDer().tau_normalized() } -> std::same_as<typename T::ScalarArray>;
        { ct.glsDer().eta_normalized() } -> std::same_as<typename T::VectorTypeArray>;
        { ct.glsDer().kappa_normalized() } -> std::same_as<typename T::ScalarArray>;
    };

    template <typename T>
    concept ProvidesGeomVar = requires(const T ct, typename T::Scalar s) {
        ct.geomVar();

        { ct.geomVar().geomVar(s, s, s) } -> std::same_as<typename T::Scalar>;
    };

    template <typename T>
    concept ProvidesPlane =
        ProvidesImplicitPrimitive<T> && requires(T t, const T ct, typename T::Scalar s, typename T::VectorType v) {
            t.plane();
            ct.plane();

            t.plane().setPlane(v, v);
        };

    template <typename T>
    concept ProvidesLine =
        ProvidesImplicitPrimitive<T> && requires(T t, const T ct, typename T::Scalar s, typename T::VectorType v) {
            t.line();
            ct.line();

            t.line().setLine(v, v);
        };

    template <typename T>
    concept ProvidesPositionCovariance = requires(const T ct) {
        ct.covarianceFit();

        { ct.covarianceFit().solver() } -> std::convertible_to<typename T::Solver>;
    };

    template <typename T>
    concept ProvidesPositionCovarianceDer = requires(const T ct) { ct.covarianceFitDer(); };

    template <typename T>
    concept ProvidesCovariancePlaneDer = requires(const T ct) {
        ct.covariancePlaneDer();

        { ct.covariancePlaneDer().dPotential() } -> std::convertible_to<typename T::ScalarArray>;
    };

    template <typename T>
    concept ProvidesNormalDer = requires(const T ct) {
        ct.covariancePlaneDer();

        { ct.covariancePlaneDer().dNormal() } -> std::convertible_to<typename T::VectorArray>;
    };

    template <typename T>
    concept ProvidesTangentPlaneBasis = requires(const T ct, typename T::VectorType v) {
        ct.tangentPlaneBasis();

        { ct.tangentPlaneBasis().worldToTangentPlane(v, true) } -> std::same_as<typename T::VectorType>;
        { ct.tangentPlaneBasis().tangentPlaneToWorld(v, true) } -> std::same_as<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesPrincipalCurvatures = requires(const T ct) {
        { ct.kmin() } -> std::same_as<typename T::Scalar>;
        { ct.kmax() } -> std::same_as<typename T::Scalar>;
        { ct.kMean() } -> std::same_as<typename T::Scalar>;
        { ct.GaussianCurvature() } -> std::same_as<typename T::Scalar>;

        { ct.kminDirection() } -> std::convertible_to<typename T::VectorType>;
        { ct.kmaxDirection() } -> std::convertible_to<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesFirstFondamentalFormComponents = requires(const T ct, typename T::Scalar s) {
        ct.firstFondamentalFormComponent();
        ct.firstFondamentalFormComponent().firstFundamentalFormComponents(s, s, s);
    };

    template <typename T>
    concept ProvidesSecondFondamentalFormComponents = requires(const T ct, typename T::Scalar s) {
        ct.secondFondamentalFormComponent();
        ct.secondFondamentalFormComponent().secondFundamentalFormComponents(s, s, s);
    };

    template <typename T>
    concept ProvidesWeingartenMap = requires(const T ct, typename T::Matrix2 m) {
        { ct.weingartenMap() } -> std::convertible_to<typename T::Matrix2>;
        ct.weingartenMap(m);
    };
} // namespace Ponca
