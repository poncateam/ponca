/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/set.h>

#include "../Common/DeviceUtils.h"
#include "ComputeObject.h"
#include "Filters.h"


/**
 * \brief Helper function to compute the output dimension of a computation
 * 
 * This function maps input shape to the corresponding outputs:
 * - (N) -> (N, _Result's shape)
 * - (N, D) -> (N, _Result's shape)
 * - (N, 1, D) -> (N, _Result's shape) (dimension collapsed)
 * - (N, K, D) -> (N, K, _Result's shape)
 * 
 * Results shape depends on the actual type. 
 * - Scalar -> () (empty, collapsed)
 * - Matrix (M, 1) -> (M,) (out dimension collapsed)
 * - Matrix (M, N) -> (M, N) (no collapsing, even if M is 1)
 * 
 * \param array The input array
 */
template<typename Point, typename _Result>
inline std::vector<size_t> ComputeVectorDim(const nb::ndarray<>& array)
{
    using Result = std::remove_cv_t<_Result>;

    if constexpr (std::is_floating_point_v<Result>)
    {
        if (array.data() == nullptr) return { 0 };
        if (array.ndim() == 1 || array.shape(1) == 1) return { array.shape(0) };
        return { array.shape(0), array.shape(1) };
    }
    else 
    {  
        static_assert(Result::ColsAtCompileTime > 0 && Result::RowsAtCompileTime > 0, "ComputeVectorDim only accepts fixed size matrices");
        
        if (array.data() == nullptr)
        {
            if constexpr (Result::ColsAtCompileTime == 1)
                return { 0, Result::RowsAtCompileTime };
            return { 0, Result::RowsAtCompileTime, Result::ColsAtCompileTime };
        }

        if (array.ndim() == 1 || array.shape(1) == 1)
        {
            if constexpr (Result::ColsAtCompileTime == 1)
                return { array.shape(0), Result::RowsAtCompileTime };
            return { array.shape(0), Result::RowsAtCompileTime, Result::ColsAtCompileTime };
        }

        // Assume this is an Eigen vector with fixed dimension (which is the case for 
        // Ponca classical return type)
        if constexpr (Result::ColsAtCompileTime == 1)
            return { array.shape(0), array.shape(1), Result::RowsAtCompileTime };
        return { array.shape(0), array.shape(1), Result::RowsAtCompileTime, Result::ColsAtCompileTime };
    }
}

/**
 * \brief Write result to the corresponding dimension
 * 
 * This is a helper function to dispatch the writting of Scalar / Vector 
 * into the output destination. 
 * 
 * \param out The output array
 * \param shift A base shift for pointer output
 * \param object The object to write
 */
template<typename Scalar, typename T>
PONCA_MULTIARCH void Write(DeviceArrayView<Scalar>& out, size_t shift, const T& object)
{
    Scalar* base = out.data + shift;

    // Note: There is no warp divergence here. As the result type of the computation will 
    // be the same for every thread. 
    if constexpr (std::is_floating_point_v<T>)
    {
        *base = object;
    }
    else if constexpr (!std::is_floating_point_v<T>) 
    {
        int strideX = (out.N > 2) ? out.stride[2] : 0;
        int strideY = (out.N > 3) ? out.stride[3] : 0;
        for (unsigned int k = 0; k < object.rows(); ++k)
            for (unsigned int l = 0; l < object.cols(); ++l)
                *(base + k * strideX + l * strideY) = object(k, l);
    }
}

/**
 * \brief Alias to a lambda function whith an input vector
 */
#define MAKE_VECTOR_FUNCTION(name) [] PONCA_MULTIARCH (CO& object, const typename CO::VectorType& in) {  return object . name ( in ); }

/**
 * \brief Helper function to run a function that takes a vector as input
 * 
 * \param object The compute object to extract result from
 * \param f The function that will be called. The signature is (object, vector)
 * \param in The input array. Will be iterated over and feed to the function
 * \param out The output location
 * \param idx Current index of computation. Allows to shift out to the correct memory location. 
 */
template<typename CO, typename Func, typename Scalar = typename CO::Scalar>
PONCA_MULTIARCH void RunVectorMethod(CO& object, Func&& f, const DeviceArrayView<Scalar>& in, DeviceArrayView<Scalar>& out, size_t idx)
{
    using Point = typename CO::DataPoint;
    // Note: there are no warp diveregence here as in is shared accross the same computation, 
    // every thread will see the same attributes.  
    
    // An array of inputs
    if (in.N == 2 && in.shape[1] == Point::Dim)
    {   
        for (size_t j = 0; j < in.shape[0]; ++j)
        {
            Write(out, idx * out.stride[0] + j * out.stride[1], f(object, in.template GetVector<Point>(j)));
        }
    }
    // A single vector
    else if (in.N ==  1 && in.shape[0] == Point::Dim)
    {
        Write(out, idx * out.stride[0], f(object, in.template AsVector<Point>()));
    }
}
