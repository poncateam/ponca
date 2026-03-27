#pragma once

#include "../Common/concepts.h"
#include "enums.h"

namespace Ponca
{
    template <typename T>
    concept ProvidesFittingDefaultTypes = ProvidesCommonTypes<T> && requires { typename T::NeighborFilter; };

    template <typename T>
    concept ProvidesPrimitiveBase = ProvidesFittingDefaultTypes<T> && requires(T t, const T ct) {
        t.setNeighborFilter(typename T::NeighborFilter{});
        { ct.getNeighborFilter() } -> std::convertible_to<typename T::NeighborFilter>;

        { ct.getCurrentState() } -> std::same_as<FIT_RESULT>;

        { ct.getWeightSum() } -> std::same_as<typename T::Scalar>;
        { ct.isStable() } -> std::same_as<bool>;
        { ct.isReady() } -> std::same_as<bool>;
        { ct.getNumNeighbors() } -> std::integral;
    };

    template <typename T>
    concept ProvidesPrimitiveDerivative = ProvidesFittingDefaultTypes<T> && requires(T t, const T ct) {
        { ct.isScaleDer() } -> std::same_as<bool>;
        { ct.isSpaceDer() } -> std::same_as<bool>;
        { ct.derDimension() } -> std::integral;
    };

    template <typename T>
    concept ProvidesAlgebraicSphere = requires(T t, const T ct, typename T::VectorType v, typename T::Scalar s) {
        ct.algebraicSphere();

        { ct.algebraicSphere().potential() } -> std::same_as<typename T::Scalar>;
        { ct.algebraicSphere().potential(v) } -> std::same_as<typename T::Scalar>;

        { ct.algebraicSphere().project(v) } -> std::same_as<typename T::VectorType>;

        { ct.algebraicSphere().primitiveGradient() } -> std::convertible_to<typename T::VectorType>;
        { ct.algebraicSphere().primitiveGradient(v) } -> std::same_as<typename T::VectorType>;

        { ct.algebraicSphere().isPlane() } -> std::same_as<bool>;
        { ct.algebraicSphere().isValid() } -> std::same_as<bool>;
        { ct.algebraicSphere().isNormalized() } -> std::same_as<bool>;
        { ct.algebraicSphere().isApprox(ct, s) } -> std::same_as<bool>;

        { t.algebraicSphere().applyPrattNorm() } -> std::same_as<bool>;
        { ct.algebraicSphere().prattNorm() } -> std::same_as<typename T::Scalar>;
        { ct.algebraicSphere().prattNorm2() } -> std::same_as<typename T::Scalar>;

        { ct.algebraicSphere().radius() } -> std::convertible_to<typename T::Scalar>;
        { ct.algebraicSphere().center() } -> std::convertible_to<typename T::VectorType>;

        t.algebraicSphere().changeBasis(v);
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
    concept ProvidesPlane = requires(T t, const T ct, typename T::Scalar s, typename T::VectorType v) {
        t.plane();
        ct.plane();

        t.plane().setPlane(v, v);
        { ct.plane().potential() } -> std::same_as<typename T::Scalar>;
        { ct.plane().potential(v) } -> std::same_as<typename T::Scalar>;

        // TODO: Add this function to plane for coherence with Algebraic sphere ?
        // { ct.plane().project() } -> std::convertible_to<typename T::VectorType>;
        { ct.plane().project(v) } -> std::same_as<typename T::VectorType>;

        { ct.plane().primitiveGradient() } -> std::same_as<typename T::VectorType>; 
        { ct.plane().primitiveGradient(v) } -> std::same_as<typename T::VectorType>; 
    };

    template<typename T>
    concept ProvidesPositionCovariance = requires(const T ct) {
        ct.covarianceFit();

        { ct.covarianceFit().solver() } -> std::convertible_to<typename T::Solver>;
    };

    template<typename T>
    concept ProvidesPositionCovarianceDer = requires(const T ct) {
        ct.covarianceFitDer();
    };


} // namespace Ponca
