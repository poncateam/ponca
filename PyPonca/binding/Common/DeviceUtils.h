/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <Ponca/Ponca>

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

namespace nb = nanobind;

#ifdef __CUDACC__
/**
 * \brief Throws an exception (visible from python) if a cuda call failed
 *
 * \param err The value returned by the cuda call
 * \param what Description of the operation that was performed
 */
inline void CheckCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess)
        throw std::runtime_error(std::string(what) + " (" + cudaGetErrorString(err) + ")");
}

/**
 * \brief Sets the current cuda device for the lifetime of the object
 *
 * Note: idea from pytorch source code https://github.com/pytorch/pytorch/blob/main/c10/cuda/CUDAGuard.h
 */
struct CudaDeviceGuard
{
    explicit CudaDeviceGuard(int device)
    {
        CheckCuda(cudaGetDevice(&m_previous), "Unable to query the current cuda device");
        CheckCuda(cudaSetDevice(device), "Unable to set the cuda device");
    }

    ~CudaDeviceGuard() { cudaSetDevice(m_previous); }

    CudaDeviceGuard(const CudaDeviceGuard&)            = delete;
    CudaDeviceGuard& operator=(const CudaDeviceGuard&) = delete;

private:
    int m_previous = 0;
};
#endif

/**
 * \brief Creates a nanobind array for return purposes
 *
 * This function is not meant to be a replacement of nb::ndarray, just a convienient way
 * to create a contiguous array.
 * The goal is to return the array to the user in the python side. It can be converted to
 * many other's library tensor with a 'lib.from_dlpack' call.
 *
 * \tparam Scalar The scalar type of the array
 *
 * \param shape Shape of the vector. Also used to describe the number of dimensions
 * \param device Device id which must match nanobind.
 * \param deviceId Index of the device (only meaningful for multi-gpu setups). Allocations
 * are performed on the current device, which must match this index.
 */
template <typename Scalar>
inline nb::ndarray<Scalar, nb::array_api> CreateDeviceArray(const std::vector<size_t>& shape, int device,
                                                            int deviceId = 0)
{
    if (shape.size() == 0)
        throw std::runtime_error("Array with shape of size 0 requested.");

    size_t elementCount = 1;
    for (unsigned int i = 0; i < shape.size(); ++i)
        elementCount *= shape[i];

    switch (device)
    {
    case nb::device::cpu::value: {
        Scalar* data = new Scalar[elementCount];
        auto deleter = nb::capsule(data, [](void* p) noexcept { delete[] reinterpret_cast<Scalar*>(p); });

        return nb::ndarray<Scalar, nb::array_api>(data, shape.size(), shape.data(), deleter, nullptr,
                                                  nb::dtype<Scalar>(), device);
    }
#ifdef __CUDACC__
    case nb::device::cuda::value: {
        // Unlike new[], cudaMalloc returns a nullptr for empty allocations, which would be
        // seen as "no output" by DeviceArrayView. We always allocate at least one element.
        Scalar* data = nullptr;
        CheckCuda(cudaMalloc(&data, std::max<size_t>(elementCount, 1) * sizeof(Scalar)), "Device allocation failed");
        auto deleter = nb::capsule(data, [](void* p) noexcept { cudaFree(p); });

        return nb::ndarray<Scalar, nb::array_api>(data, shape.size(), shape.data(), deleter, nullptr,
                                                  nb::dtype<Scalar>(), device, deviceId);
    }
#endif
    default:
        throw std::runtime_error("Unsupported device for allocation");
    }
}

/**
 * \brief Creates an array on a device
 *
 * Note: this is a bitwise copy. As nanobind can only hold arithmetic types, the results
 * is within bytes to be reinterpreted as T*.
 *
 * \param elements The host elements
 * \param device Device id which must match nanobind.
 * \param deviceId Index of the device (see CreateDeviceArray)
 */
template <typename T>
inline nb::ndarray<uint8_t, nb::array_api> CopyToDevice(const std::vector<T>& elements, int device, int deviceId = 0)
{
    const size_t bytes = elements.size() * sizeof(T);
    auto array         = CreateDeviceArray<uint8_t>({bytes}, device, deviceId);
    if (bytes == 0)
        return array;

    switch (device)
    {
    case nb::device::cpu::value:
        std::memcpy(array.data(), elements.data(), bytes);
        break;
#ifdef __CUDACC__
    case nb::device::cuda::value:
        CheckCuda(cudaMemcpy(array.data(), elements.data(), bytes, cudaMemcpyHostToDevice), "Copy to device failed");
        break;
#endif
    }
    return array;
}

#ifdef __CUDACC__
/**
 * \brief Copies a device array to the host
 *
 * \param array The device array
 */
template <typename Scalar, typename... Args>
nb::ndarray<Scalar, Args...> CopyToHost(const nb::ndarray<Scalar, Args...>& array)
{
    const auto elementStride = [](int64_t stride) {
        return static_cast<int64_t>(static_cast<uint64_t>(stride) * sizeof(Scalar)) /
               static_cast<int64_t>(sizeof(Scalar));
    };

    int64_t first = 0, last = 0;
    std::vector<size_t> shape(array.ndim());
    std::vector<int64_t> strides(array.ndim());
    for (size_t d = 0; d < array.ndim(); ++d)
    {
        shape[d]   = array.shape(d);
        strides[d] = elementStride(array.stride(d));

        const int64_t extent = (static_cast<int64_t>(array.shape(d)) - 1) * strides[d];
        first += std::min<int64_t>(0, extent);
        last += std::max<int64_t>(0, extent);
    }

    const size_t count = (array.size() == 0) ? 0 : static_cast<size_t>(last - first + 1);
    Scalar* span       = new Scalar[std::max<size_t>(count, 1)];
    auto deleter       = nb::capsule(span, [](void* p) noexcept { delete[] reinterpret_cast<Scalar*>(p); });

    CheckCuda(cudaMemcpy(span, array.data() + first, count * sizeof(Scalar), cudaMemcpyDeviceToHost),
              "Copy to host failed");

    return nb::ndarray<Scalar, Args...>(span - first, shape.size(), shape.data(), deleter, strides.data(),
                                        nb::dtype<Scalar>(), nb::device::cpu::value);
}
#endif

/**
 * \brief Array view over a nanobind array accessible on device
 *
 * This class creates a view of an array, pushing strides and shapes
 * on the target device. It DOESN't TAKE OWNERSHIP OF DATA NOR DOES IT
 * COPIES IT. It only copies shape and stride to the destination device.
 *
 * It serves as a convenient way to always have the same data structure
 * no matter the backend.
 *
 * This class is not meant as a replacement of nb::ndarray, especially since
 * shape and strides may not be accessible in host code.
 */
template <typename Scalar>
struct DeviceArrayView
{
    /**
     * \brief Default ctor
     */
    DeviceArrayView() : shape(nullptr), stride(nullptr)
    {
        N      = 0;
        data   = nullptr;
        device = -1;
    }

    /**
     * \brief Constructor from a nb array
     *
     * Note: shape and stride will only be available on the device
     * on the same device as data
     */
    DeviceArrayView(nanobind::ndarray<Scalar> array) : shape(nullptr), stride(nullptr)
    {
        data   = array.data();
        N      = array.ndim();
        device = array.device_type();

        Copy(nb::device::cpu::value, array.shape_ptr(), array.stride_ptr());
    }

    /**
     * \brief Copy constructor
     *
     * A copy of shape and stride is performed, the ownership is not shared.
     * This may (or may not) support device conversion.
     */
    DeviceArrayView(const DeviceArrayView& other) : shape(nullptr), stride(nullptr)
    {
        data   = other.data;
        N      = other.N;
        device = other.device;

        Copy(other.device, other.shape, other.stride);
    }

    /**
     * \brief Assignment
     */
    DeviceArrayView& operator=(const DeviceArrayView& other)
    {
        if (&other != this)
        {
            data   = other.data;
            N      = other.N;
            device = other.device;

            Copy(other.device, other.shape, other.stride);
        }
        return *this;
    }

    /**
     * \brief Convert the array into a Vector
     *
     * This function is mean for shape like (D, ) i.e. the array describes
     * a single vector. It does not check if data, stride or anything else is valid.
     * It just iterates the first available dimension, filling a vector as
     * much as possible.
     *
     * \see GetVector
     */
    template <typename Point>
    PONCA_MULTIARCH typename Point::VectorType AsVector() const
    {
        typename Point::VectorType out;
        for (size_t j = 0; j < shape[0]; ++j)
            out(j) = *(data + j * stride[0]);
        return out;
    }

    /**
     * \brief Fill a given array
     *
     * This function is mean for shape like (D, ) i.e. the array describes
     * a single vector. It does not check if data, stride or anything else is valid.
     * It just iterates the first available dimension, filling a vector as
     * much as possible.
     *
     * \see GetVector
     */
    PONCA_MULTIARCH void Fill(Scalar* out) const
    {
        for (size_t j = 0; j < shape[0]; ++j)
            out[j] = *(data + j * stride[0]);
    }

    /**
     * \brief Convert the array into a Vector
     *
     * This function is meant for shape like (N, D) i.e. the array describes
     * an array of vectors. It does not check if data, stride or anything else
     * is valid. It just iterates the first available dimension, filling a vector as
     * much as possible.
     *
     * As a side note, this function can not be used to convert an array of shape (D, )
     * to a Vector as the function assumes that the shape as two components...
     *
     * \see GetVector
     */
    template <typename Point>
    PONCA_MULTIARCH typename Point::VectorType GetVector(size_t i) const
    {
        typename Point::VectorType out;
        for (size_t j = 0; j < shape[1]; ++j)
            out(j) = *(data + i * stride[0] + j * stride[1]);
        return out;
    }

    /**
     * \brief Fills an an array
     *
     * This function is meant for shape like (N, D) i.e. the array describes
     * an array of vectors. It does not check if data, stride or anything else
     * is valid. It just iterates the first available dimension, filling a vector as
     * much as possible.
     *
     * As a side note, this function can not be used to convert an array of shape (D, )
     * to a Vector as the function assumes that the shape as two components...
     *
     * \see GetVector
     */
    PONCA_MULTIARCH void Fill(Scalar* out, size_t i) const
    {
        for (size_t j = 0; j < shape[1]; ++j)
            out[j] = *(data + i * stride[0] + j * stride[1]);
    }

    /**
     * \brief Convert the array into a Vector
     *
     * This function is meant for shape like (N, K, D) i.e. the array describes
     * an array of arrays of vectors. It does not check if data, stride or anything else
     * is valid. It just iterates the dimensions, filling a vector as much as possible.
     *
     * \see GetVector
     */
    template <typename Point>
    PONCA_MULTIARCH typename Point::VectorType GetVector(size_t i, size_t j) const
    {
        typename Point::VectorType out;
        for (size_t k = 0; k < shape[2]; ++k)
            out(k) = *(data + i * stride[0] + j * stride[1] + k * stride[2]);
        return out;
    }

    /**
     * \brief Fills an an array
     *
     * This function is meant for shape like (N, D) i.e. the array describes
     * an array of vectors. It does not check if data, stride or anything else
     * is valid. It just iterates the first available dimension, filling a vector as
     * much as possible.
     *
     * As a side note, this function can not be used to convert an array of shape (D, )
     * to a Vector as the function assumes that the shape as two components...
     *
     * \see GetVector
     */
    PONCA_MULTIARCH void Fill(Scalar* out, size_t i, size_t j) const
    {
        for (size_t k = 0; k < shape[2]; ++k)
            out[k] = *(data + i * stride[0] + j * stride[1] + k * stride[2]);
    }

    /**
     * \brief Free data
     */
    ~DeviceArrayView() { Free(); }

    /**
     * \brief Return device
     *
     * Note: This is function is meant to be compatible with generic code mixing
     * nb::ndarray and this.
     */
    int device_type() const { return device; }

    // We leave free access to data, this is just a copy, if the user wants to mess
    // with it, then so be it.

    size_t N;
    int device;
    Scalar* data;

    int64_t* shape;
    int64_t* stride;

private:
    /**
     * \brief Free shape and stride according to the device they are on
     */
    void Free()
    {
        if (shape == nullptr || stride == nullptr)
            return;

        switch (device)
        {
        case nb::device::cpu::value: {
            delete[] shape;
            delete[] stride;
        }
        break;
#ifdef __CUDACC__
        case nb::device::cuda::value: {
            cudaFree(shape);
            cudaFree(stride);
        }
        break;
#endif
        default:
            throw std::runtime_error("Unsupported device when freeing (" + std::to_string(device) + ")!");
        }

        shape  = nullptr;
        stride = nullptr;
    }

    /**
     * \brief Copy shape and stride to the corresponding device
     *
     * \param srcDevice The device srcshape and srcstride are located on
     * \param srcshape A pointer to the shape data
     * \param srcstride A pointer to the stride data
     */
    void Copy(int srcDevice, const int64_t* srcshape, const int64_t* srcstride)
    {
        // Note: We could pack shape and stride into a single memory calls
        // but for now this is simpler :)
        Free();
        if (data == nullptr)
            return;

        if (srcshape == nullptr && srcstride == nullptr)
            return;

        switch (device)
        {
        case nb::device::cpu::value: {
            shape  = new int64_t[N];
            stride = new int64_t[N];

            switch (srcDevice)
            {
            case nb::device::cpu::value:
                memcpy(shape, srcshape, N * sizeof(int64_t));
                memcpy(stride, srcstride, N * sizeof(int64_t));
                break;
#ifdef __CUDACC__
            case nb::device::cuda::value:
                cudaMemcpy(shape, srcshape, sizeof(size_t) * N, cudaMemcpyDeviceToHost);
                cudaMemcpy(stride, srcstride, sizeof(size_t) * N, cudaMemcpyDeviceToHost);
                break;
#endif
            default:
                Free();
                throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" +
                                         std::to_string(srcDevice) + ")!");
            }
        }
        break;
#ifdef __CUDACC__
        case nb::device::cuda::value: {
            cudaMalloc(&shape, sizeof(int64_t) * N);
            cudaMalloc(&stride, sizeof(int64_t) * N);

            cudaMemcpyKind cpyMode;
            switch (srcDevice)
            {
            case nb::device::cpu::value:
                cpyMode = cudaMemcpyHostToDevice;
                break;
            case nb::device::cuda::value:
                cpyMode = cudaMemcpyDeviceToDevice;
                break;
            default:
                Free();
                throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" +
                                         std::to_string(srcDevice) + ")!");
            }

            cudaMemcpy(shape, srcshape, sizeof(int64_t) * N, cpyMode);
            cudaMemcpy(stride, srcstride, sizeof(int64_t) * N, cpyMode);
        }
        break;
#endif
        default:
            Free();
            throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" +
                                     std::to_string(srcDevice) + ")!");
        }
    }
};
