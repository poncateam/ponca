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
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <stdexcept>
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
 * \brief Sets the current device for the lifetime of the object, whatever the device
 *
 * Does nothing for the CPU, and when compiled without CUDA support. Allocations (see AllocateBytes)
 * already select their device: this guard is only needed to launch kernels.
 *
 * \param device Device type, which must match nanobind.
 * \param deviceId Index of the device
 */
struct DeviceGuard
{
    DeviceGuard(int device, int deviceId)
    {
#ifdef __CUDACC__
        if (device == nb::device::cuda::value)
            m_guard.emplace(deviceId);
#else
        (void)device;
        (void)deviceId;
#endif
    }

    DeviceGuard(const DeviceGuard&)            = delete;
    DeviceGuard& operator=(const DeviceGuard&) = delete;

#ifdef __CUDACC__
private:
    std::optional<CudaDeviceGuard> m_guard;
#endif
};

/**
 * \brief Throws if the given device is not the CPU, for operations that are only available on the host
 *
 * \param device Device type, which must match nanobind.
 * \param what Description of the operation, used as error message
 */
inline void RequireHost(int device, const char* what)
{
    if (device != nb::device::cpu::value)
        throw std::runtime_error(std::string(what) + " (got device type " + std::to_string(device) + ")");
}

/**
 * \brief Allocates memory on a device
 *
 * At least one byte is always allocated, so that a nullptr can be seen as "no data".
 * The memory must be released with FreeBytes.
 *
 * \param bytes Size of the allocation
 * \param device Device type, which must match nanobind.
 * \param deviceId Index of the device
 */
inline void* AllocateBytes(size_t bytes, int device, int deviceId = 0)
{
    bytes = std::max<size_t>(bytes, 1);
    switch (device)
    {
    case nb::device::cpu::value:
        return ::operator new(bytes);
#ifdef __CUDACC__
    case nb::device::cuda::value: {
        DeviceGuard guard(device, deviceId);
        void* data = nullptr;
        CheckCuda(cudaMalloc(&data, bytes), "Device allocation failed");
        return data;
    }
#endif
    default:
        throw std::runtime_error("Unsupported device for allocation (" + std::to_string(device) + ")");
    }
}

/**
 * \brief Releases memory allocated by AllocateBytes
 *
 * \param data The memory, which may be nullptr
 * \param device Device type the memory was allocated on
 */
inline void FreeBytes(void* data, int device) noexcept
{
    if (data == nullptr)
        return;

    switch (device)
    {
    case nb::device::cpu::value:
        ::operator delete(data);
        break;
#ifdef __CUDACC__
    case nb::device::cuda::value:
        cudaFree(data);
        break;
#endif
    default:
        break; // AllocateBytes does not support other devices
    }
}

/**
 * \brief Creates a capsule releasing memory allocated by AllocateBytes, to own nanobind arrays
 *
 * \param data The memory
 * \param device Device type the memory was allocated on
 */
inline nb::capsule MakeOwner(void* data, int device)
{
    // Capsule deleters can not capture the device
    constexpr auto freeCpu = [](void* p) noexcept { FreeBytes(p, nb::device::cpu::value); };
#ifdef __CUDACC__
    constexpr auto freeCuda = [](void* p) noexcept { FreeBytes(p, nb::device::cuda::value); };
#endif

    switch (device)
    {
    case nb::device::cpu::value:
        return nb::capsule(data, freeCpu);
#ifdef __CUDACC__
    case nb::device::cuda::value:
        return nb::capsule(data, freeCuda);
#endif
    default:
        FreeBytes(data, device);
        throw std::runtime_error("Unsupported device for allocation (" + std::to_string(device) + ")");
    }
}

/**
 * \brief Copies memory between two devices (in any direction)
 *
 * \param dst Destination
 * \param dstDevice Device type of the destination
 * \param src Source
 * \param srcDevice Device type of the source
 * \param bytes Number of bytes to copy
 */
inline void CopyBytes(void* dst, int dstDevice, const void* src, int srcDevice, size_t bytes)
{
    if (bytes == 0)
        return;

    constexpr int cpu = nb::device::cpu::value;
    if (srcDevice == cpu && dstDevice == cpu)
    {
        std::memcpy(dst, src, bytes);
        return;
    }

#ifdef __CUDACC__
    constexpr int cuda = nb::device::cuda::value;
    if ((srcDevice == cpu || srcDevice == cuda) && (dstDevice == cpu || dstDevice == cuda))
    {
        const cudaMemcpyKind kind = srcDevice == cpu   ? cudaMemcpyHostToDevice
                                    : dstDevice == cpu ? cudaMemcpyDeviceToHost
                                                       : cudaMemcpyDeviceToDevice;
        CheckCuda(cudaMemcpy(dst, src, bytes, kind), "Copy between devices failed");
        return;
    }
#endif

    throw std::runtime_error("Unsupported copy from device " + std::to_string(srcDevice) + " to device " +
                             std::to_string(dstDevice));
}

/**
 * \brief Computes the memory span covered by a (strided) array
 *
 * Strides are given in elements, and may be 0 (broadcasted dimensions) or negative.
 *
 * \param array The array
 * \param first Offset of the first element of the span from the array data (<= 0)
 * \param count Number of elements in the span
 */
template <typename Array>
void ComputeArraySpan(const Array& array, int64_t& first, size_t& count)
{
    first        = 0;
    int64_t last = 0;
    for (size_t d = 0; d < array.ndim(); ++d)
    {
        const int64_t extent = (static_cast<int64_t>(array.shape(d)) - 1) * array.stride(d);
        first += std::min<int64_t>(0, extent);
        last += std::max<int64_t>(0, extent);
    }
    count = (array.size() == 0) ? 0 : static_cast<size_t>(last - first + 1);
}

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
 * \param device Device type, which must match nanobind.
 * \param deviceId Index of the device (only meaningful for multi-gpu setups).
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

    Scalar* data = static_cast<Scalar*>(AllocateBytes(elementCount * sizeof(Scalar), device, deviceId));
    return nb::ndarray<Scalar, nb::array_api>(data, shape.size(), shape.data(), MakeOwner(data, device), nullptr,
                                              nb::dtype<Scalar>(), device, deviceId);
}

/**
 * \brief Copies an array to a device (in any direction), keeping its shape and strides
 *
 * Only the memory span covered by the array is copied, so broadcasted (0-stride)
 * dimensions are not expanded. The returned array owns its memory.
 *
 * \param array The array to copy
 * \param device Device type, which must match nanobind.
 * \param deviceId Index of the device
 */
template <typename Scalar, typename... Args>
nb::ndarray<Scalar, Args...> CopyArray(const nb::ndarray<Scalar, Args...>& array, int device, int deviceId = 0)
{
    int64_t first;
    size_t count;
    ComputeArraySpan(array, first, count);
    const size_t bytes = count * sizeof(Scalar);

    Scalar* storage   = static_cast<Scalar*>(AllocateBytes(bytes, device, deviceId));
    nb::capsule owner = MakeOwner(storage, device);
    CopyBytes(storage, device, array.data() + first, array.device_type(), bytes);

    const std::vector<size_t> shape(array.shape_ptr(), array.shape_ptr() + array.ndim());
    return nb::ndarray<Scalar, Args...>(storage - first, shape.size(), shape.data(), owner, array.stride_ptr(),
                                        nb::dtype<Scalar>(), device, deviceId);
}

/**
 * \brief Typed, reference counted buffer of elements stored on a device
 *
 * Elements are copied bitwise: T must be valid on the target device once copied (e.g. pointers
 * it holds must point to memory of that device). Copies of the buffer share the same memory.
 *
 * \tparam T Type of the elements
 */
template <typename T>
class DeviceBuffer
{
public:
    DeviceBuffer() = default;

    /**
     * \brief Copies elements stored on any device
     *
     * \param elements The elements to copy
     * \param size Number of elements
     * \param srcDevice Device type the elements are stored on
     * \param device Device type of the buffer, which must match nanobind.
     * \param deviceId Index of the device
     */
    DeviceBuffer(const T* elements, size_t size, int srcDevice, int device, int deviceId = 0)
        : m_size(size), m_device(device), m_deviceId(deviceId)
    {
        const size_t bytes = size * sizeof(T);
        void* data         = AllocateBytes(bytes, device, deviceId);
        m_storage          = std::shared_ptr<void>(data, [device](void* p) noexcept { FreeBytes(p, device); });
        CopyBytes(data, device, elements, srcDevice, bytes);
    }

    /**
     * \brief Copies elements stored on the host
     *
     * \param elements The elements to copy
     * \param device Device type of the buffer, which must match nanobind.
     * \param deviceId Index of the device
     */
    DeviceBuffer(const std::vector<T>& elements, int device, int deviceId = 0)
        : DeviceBuffer(elements.data(), elements.size(), nb::device::cpu::value, device, deviceId)
    {
    }

    /**
     * \brief Returns a copy of the buffer, stored on another device
     */
    DeviceBuffer toDevice(int device, int deviceId = 0) const
    {
        return DeviceBuffer(data(), m_size, m_device, device, deviceId);
    }

    T* data() const { return static_cast<T*>(m_storage.get()); }

    size_t size() const { return m_size; }

    int device_type() const { return m_device; }

    int device_id() const { return m_deviceId; }

private:
    std::shared_ptr<void> m_storage;
    size_t m_size{0};
    int m_device{nb::device::cpu::value};
    int m_deviceId{0};
};

/**
 * \brief Array view over a nanobind array accessible on device
 *
 * This class creates a view of an array, storing its shape and strides inline.
 * It DOESN'T TAKE OWNERSHIP OF DATA NOR DOES IT COPY IT.
 *
 * It is trivially copyable, so it can be copied bytewise to any device (see DeviceBuffer)
 * and serves as a convenient way to always have the same data structure no matter the backend.
 *
 * This class is not meant as a replacement of nb::ndarray.
 */
template <typename Scalar>
struct DeviceArrayView
{
    /// Maximum number of dimensions of the viewed arrays
    static constexpr size_t MaxDims = 4;

    /**
     * \brief Default ctor, viewing no data
     */
    DeviceArrayView() = default;

    /**
     * \brief Constructor from a nb array
     *
     * \throw std::runtime_error if the array has more than MaxDims dimensions
     */
    explicit DeviceArrayView(const nanobind::ndarray<Scalar>& array)
        : N(array.ndim()), device(array.device_type()), deviceId(array.device_id()), data(array.data())
    {
        if (N > MaxDims)
            throw std::runtime_error("Arrays with more than " + std::to_string(MaxDims) +
                                     " dimensions are not supported (got " + std::to_string(N) + ")");

        for (size_t d = 0; d < N; ++d)
        {
            shape[d]  = array.shape(d);
            stride[d] = array.stride(d);
        }
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
     * \brief Return device
     *
     * Note: This is function is meant to be compatible with generic code mixing
     * nb::ndarray and this.
     */
    int device_type() const { return device; }

    // We leave free access to data, this is just a copy, if the user wants to mess
    // with it, then so be it.

    size_t N     = 0;
    int device   = -1;
    int deviceId = 0;
    Scalar* data = nullptr;

    int64_t shape[MaxDims]  = {};
    int64_t stride[MaxDims] = {};
};
