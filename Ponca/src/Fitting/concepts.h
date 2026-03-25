#pragma once

#include "../Common/concepts.h"
#include "enums.h"

namespace Ponca
{
    template <typename T>
    concept ProvidesFittingDefaultTypes = ProvidesCommonTypes<T> && requires {
        typename T::NeighborFilter;
    };

    template <typename T>
    concept ProvidesPrimitiveBase =
        ProvidesFittingDefaultTypes<T> && requires(T t, const T ct) {
            t.setNeighborFilter(typename T::NeighborFilter{});
            { ct.getNeighborFilter() } -> std::convertible_to<typename T::NeighborFilter>;

            { ct.getCurrentState() } -> std::same_as<FIT_RESULT>;
            
            { ct.getWeightSum() } -> std::same_as<typename T::Scalar>;
            { ct.isStable() } -> std::same_as<bool>;
            { ct.isReady() } -> std::same_as<bool>;
            { ct.getNumNeighbors() } -> std::integral;
        };

    template <typename T>
    concept ProvidesPrimitiveDerivative = 
        ProvidesFittingDefaultTypes<T> && requires(T t, const T ct) {
            { ct.isScaleDer() } -> std::same_as<bool>;
            { ct.isSpaceDer() } -> std::same_as<bool>;
            { ct.derDimension() } -> std::integral;
        };

} // namespace Ponca
