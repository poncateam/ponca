/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <Ponca/Ponca>
#include <type_traits>

namespace nb = nanobind;

/**
 * \brief Equivalent to a list of Scalar
 */
template <typename Point>
using PyScalarArray = nb::ndarray<typename Point::Scalar, nb::shape<-1>>;

/**
 * \brief Equivalent to typename Point::VectorType
 */
template <typename Point>
using PyVector = nb::ndarray<typename Point::Scalar, nb::shape<Point::Dim>>;

/**
 * \brief Equivalent to an array of typename Point::VectorType
 */
template <typename Point>
using PyVectorArray = nb::ndarray<typename Point::Scalar, nb::shape<-1, Point::Dim>>;

/**
 * \brief Equivalent to an array of arrays of typename Point::VectorType
 */
template <typename Point>
using PyVectorVectorArray = nb::ndarray<typename Point::Scalar, nb::shape<-1, -1, Point::Dim>>;

/**
 * \brief Convert a Vector from python to typename Point::VectorType
 *
 * \tparam Point Point type
 *
 * \param vector A python view over a vector
 */
template <typename Point>
inline typename Point::VectorType PyVectorToVector(const PyVector<Point>& vector)
{
    using Map = Eigen::Map<typename Point::VectorType>;
    return Map(vector.data(), Point::Dim, 1);
}

/**
 * \brief Index into a Scalar array
 *
 * This function does not perform bound checking
 *
 * \tparam Point Point type
 *
 * \param array A python view over the array
 * \param idx The index
 */
template <typename Point>
inline typename Point::Scalar PyScalarArrayIndex(const PyScalarArray<Point>& array, unsigned int idx)
{
    return *(array.data() + idx * array.stride(0));
}

/**
 * \brief Index into an array of vectors
 *
 * This function does not perform bound checking
 *
 * \tparam Point Point type
 *
 * \param array A python view over the array
 * \param idx The index
 */
template <typename Point>
inline auto PyVectorArrayIndex(const PyVectorArray<Point>& array, unsigned int idx)
{
    using Map = Eigen::Map<typename Point::VectorType>;
    return Map(array.data() + idx * array.stride(0), Point::Dim, 1);
}

/**
 * \brief Index into an array of arrays of vectors
 *
 * This function does not perform bound checking
 *
 * \tparam Point Point type
 *
 * \param array A python view over the array
 * \param i First index
 * \param j Second index
 */
template <typename Point>
inline auto PyVectorVectorArrayIndex(const PyVectorVectorArray<Point>& array, unsigned int i, unsigned int j)
{
    using Map = Eigen::Map<typename Point::VectorType>;
    return Map(array.data() + i * array.stride(0) + j * array.stride(1), Point::Dim, 1);
}

/**
 * \brief Shortcut to nb capsule for standard deletter
 */
template <typename T>
inline nb::capsule PyStandardDeleter(T* p)
{
    return nb::capsule(p, [](void* p) noexcept { delete reinterpret_cast<T*>(p); });
}

/**
 * \brief Shortcut to nb capsule for standard deletter
 */
template <typename T>
inline nb::capsule PyArrayDeleter(T* p)
{
    return nb::capsule(p, [](void* p) noexcept { delete[] reinterpret_cast<T*>(p); });
}

/**
 * \brief Broadcast an array to a desired shape
 *
 * This function does not check if the broadcast can be valid !
 * THe returned array shares the data and no copies are performed.
 *
 * \tparam T Scalar type (not used, here for compatibility)
 * \tparam N Number of dimensions
 * \tparam Array Array type
 *
 * \param array The array to boradcast
 * \param targetShape The target shape
 */
template <typename T, size_t N, typename Array>
auto broadcastArrayTo(const Array& array, const std::array<size_t, N>& targetShape)
{
    const size_t ndim = array.ndim();

    if (ndim > N)
        throw std::runtime_error("Cannot broadcast to fewer dimensions");

    std::array<std::int64_t, N> strides{};
    for (size_t i = 0; i < ndim; ++i)
    {
        size_t src = ndim - 1 - i;
        size_t dst = N - 1 - i;

        auto srcExtent = array.shape(src);
        auto dstExtent = targetShape[dst];

        strides[dst] = (srcExtent == dstExtent) * array.stride(src);
    }

    return nb::ndarray<T>(array.data(), N, targetShape.data(),
                          const_cast<Array&>(array).cast(), // .cast() is not a const method...
                          strides.data(), array.dtype(), array.device_type(), array.device_id());
}
/**
 * \brief Broadcast a scalar to an array of the desired shape
 *
 * This function allocates a new single scalar.
 *
 * \tparam T Scalar type (not used, here for compatibility)
 * \tparam N Number of dimensions
 *
 * \param t The scalar to boradcast
 * \param targetShape The target shape
 */
template <typename T, size_t N>
auto broadcastTo(const T& t, const std::array<size_t, N>& targetShape)
{
    T* newdata = new T(t);
    return nb::ndarray<T>(newdata, N, targetShape.data(), PyStandardDeleter(newdata),
                          std::array<std::int64_t, N>{}.data(), nb::dtype<T>(), nb::device::cpu());
}

