#pragma once

#include "../concepts.h"

namespace Ponca
{
    /// \brief This concept if currently only used a base concept for ProvidesBasketUnitBase
    ///
    /// \fixme Find a way to make the methods getCurrentState, isStable and isReady accessible outside of the Basket.
    ///        This could be done by moving their definition to ComputeObject, unfortunately there is a design flaw
    ///        that prevent to do it currently: both the Basket and BasketDiff inherit from ComputeObject, and so any
    ///        virtual function (e.g. getCurrentState is virtual and the other are in ComputeObject only) are found
    ///        in multiple base classes of different types by the compiler.
    /// \fixme For the same reason as abobe, this concept cannot be tested out of the Basket hierarchy (e.g. as a
    ///        requirement of MLSEvaluationScheme::computeMLSImpl (which would make sense as we need access to
    ///        isStable in the MLS iterations)
    template <typename T>
    concept ProvidesComputeState = ProvidesCommonTypes<T> && requires(T t, const T ct) {
        { ct.isStable() } -> std::same_as<bool>;
        { ct.isReady() } -> std::same_as<bool>;
        { ct.getCurrentState() } -> std::same_as<FIT_RESULT>;
    };
    /// \brief This concept ensures that the default types and accessors in a Basket are well-formed
    ///
    /// This concept is implemented by BasketUnitBase, which is set by default in the Basket.
    template <typename T>
    concept ProvidesBasketUnitBase = ProvidesCommonTypes<T> && ProvidesComputeState<T> && requires(T t, const T ct) {
        typename T::NeighborFilter;
        t.setNeighborFilter(typename T::NeighborFilter{});
        { ct.getNeighborFilter() } -> std::convertible_to<typename T::NeighborFilter>;

        { ct.getWeightSum() } -> std::same_as<typename T::Scalar>;
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

} // namespace Ponca
