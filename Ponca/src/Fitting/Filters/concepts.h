#pragma once

#include "../concepts.h"

namespace Ponca
{
    template <typename T>
    concept HasLocalFrame = T::hasLocalFrame;

    template <typename T>
    concept IsNeighborhoodFrame = requires(T t, const T ct, typename T::VectorType v) {
        t.changeNeighborhoodFrame(v);

        { ct.convertToGlobalBasis(v, true) } -> std::convertible_to<typename T::VectorType>;
        { ct.convertToLocalBasis(v, true) } -> std::convertible_to<typename T::VectorType>;

        { t.center() } -> std::convertible_to<typename T::VectorType>;
        { ct.center() } -> std::convertible_to<typename T::VectorType>;
    };

    template <typename T>
    concept IsNeighborhoodFilter = requires(T t, const T ct, typename T::VectorType v, typename T::DataPoint d) {
        t.frame();
        ct.frame();

        ct.operator()(d); // not sure how to test the return type (std::pair) without breaking Cuda compatibility
        { ct.spacedw(v, d) } -> std::convertible_to<typename T::VectorType>;
        { ct.spaced2w(v, d) } -> std::convertible_to<typename T::MatrixType>;
        { ct.scaledw(v, d) } -> std::same_as<typename T::Scalar>;
        { ct.scaled2w(v, d) } -> std::same_as<typename T::Scalar>;
        { ct.scaleSpaced2w(v, d) } -> std::convertible_to<typename T::VectorType>;
    };

    template <typename T>
    concept ProvidesNeighborhoodFilter = ProvidesCommonTypes<T> && IsNeighborhoodFilter<typename T::NeighborFilter> &&
                                         IsNeighborhoodFrame<typename T::NeighborFilter> && requires(T t, const T ct) {
                                             t.getNeighborFilter();
                                             ct.getNeighborFilter();
                                         };
} // namespace Ponca
