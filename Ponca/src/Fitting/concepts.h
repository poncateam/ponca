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
} // namespace Ponca
