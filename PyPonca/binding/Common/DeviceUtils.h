/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <Ponca/Ponca>

namespace nb = nanobind;
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
 */
template <typename Scalar>
inline nb::ndarray<Scalar, nb::array_api> CreateDeviceArray(const std::vector<size_t>& shape, int device)
{
    std::vector<int64_t> strides(shape.size());
    
    strides[strides.size() - 1] = 1;
    for (int64_t i = strides.size() - 2; i >= 0; --i)
        strides[i] = strides[i - 1] * shape[i - 1];

    const size_t elementCount = strides[0] * shape[0];
    
    switch (device)
    {
        case nb::device::cpu::value:
        {
            Scalar* data = new Scalar[elementCount];
            auto deleter = nb::capsule(data, [](void* p) noexcept { delete[] reinterpret_cast<Scalar*>(p); });

            return nb::ndarray<Scalar, nb::array_api>(
                data, shape.size(), shape.data(), deleter, strides.data(), nb::dtype<Scalar>(), device
            );
        }
#ifdef __CUDACC__
        case nb::device::cuda::value:
        {
            Scalar* data = nullptr;
            cudaMalloc(&data, elementCount * sizeof(data));
            auto deleter = nb::capsule(data, [](void* p) noexcept { cudaFree(p); });

            return nb::ndarray<Scalar, nb::array_api>(
                data, shape.size(), shape.data(), deleter, strides.data(), nb::dtype<Scalar>(), device
            );
        }
#endif
        default:
            throw std::runtime_error("Unsupported device for allocation");
    }
}

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
 * 
 * TODO: Group allocations (at least of shape / stride) into a single ptr? 
 */
template<typename Scalar>
struct DeviceArrayView
{
    /**
     * \brief Default ctor
     */
    DeviceArrayView() :
        shape(nullptr), stride(nullptr)
    {
        N = 0;
        data = nullptr;
        device = -1;
    }

    /**
     * \brief Constructor from a nb array
     * 
     * Note: shape and stride will only be available on the device
     * on the same device as data 
     */
    DeviceArrayView(nanobind::ndarray<Scalar> array) : 
        shape(nullptr), stride(nullptr)
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
    DeviceArrayView(const DeviceArrayView& other) : 
        shape(nullptr), stride(nullptr)
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
    template<typename Point>
    PONCA_MULTIARCH typename Point::VectorType AsVector() const
    {
        typename Point::VectorType out;
        for (size_t j = 0; j < shape[0]; ++j)
            out(j) = *(data + j * stride[0]);
        return out;    
    }

    /**
     * \brief Convert the array into a Vector
     * 
     * This function is mean for shape like (N, D) i.e. the array describes 
     * an array of vectors. It does not check if data, stride or anything else 
     * is valid. It just iterates the first available dimension, filling a vector as 
     * much as possible.
     * 
     * As a side note, this function can not be used to convert an array of shape (D, )
     * to a Vector as the function assumes that the shape as two components... 
     * 
     * \see GetVector 
     */
    template<typename Point>
    PONCA_MULTIARCH typename Point::VectorType GetVector(size_t i) const
    {
        typename Point::VectorType out;
        for (size_t j = 0; j < shape[1]; ++j)
            out(j) = *(data + i * stride[0] + j * stride[1]);
        return out;
    }

    /**
     * \brief Free data
     */
    ~DeviceArrayView()
    {
        Free();
    }

    /**
     * \brief Return device 
     * 
     * Note: This is function is meant to be compatible with generic code mixing 
     * nb::ndarray and this. 
     */
    int device_type() const
    {
        return device;
    }

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
        if (data == nullptr) return;

        switch (device)
        {
        case nb::device::cpu::value:
            {
                delete[] shape;
                delete[] stride;
            }
            break;
#ifdef __CUDACC__
        case nb::device::cuda::value:
            {
                cudaFree(shape);
                cudaFree(stride);
            }
            break;
#endif 
        default:
            throw std::runtime_error("Unsupported device when freeing (" + std::to_string(device) + ")!");
        }

        shape = nullptr;
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

        if (srcshape == nullptr && srcstride == nullptr)
            return;

        switch(device)
        {
        case nb::device::cpu::value:
            {
                shape  = new int64_t[N];
                stride = new int64_t[N];

                switch (srcDevice)
                {
                case nb::device::cpu::value:
                    memcpy(shape , srcshape , N * sizeof(int64_t));
                    memcpy(stride, srcstride, N * sizeof(int64_t));
                    break;                    
#ifdef __CUDACC__
                case nb::device::cuda::value:
                    cudaMemcpy(shape , srcshape, sizeof(size_t) * N, cudaMemcpyDeviceToHost);
                    cudaMemcpy(stride, srcstride, sizeof(size_t) * N, cudaMemcpyDeviceToHost);
                    break;
#endif
                default:
                    Free();
                    throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" + std::to_string(srcDevice) + ")!");
                }
            }
            break;
#ifdef __CUDACC__
        case nb::device::cuda::value:
            {
                cudaMalloc(&shape , sizeof(int64_t) * N);
                cudaMalloc(&stride, sizeof(int64_t) * N);

                cudaMemcpyKind cpyMode;
                switch(srcDevice)
                {
                case nb::device::cpu::value:
                    cpyMode = cudaMemcpyHostToDevice;
                    break;
                case nb::device::cuda::value:
                    cpyMode = cudaMemcpyHostToHost;
                    break;
                default:
                    Free();
                    throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" + std::to_string(srcDevice) + ")!");
                }
                
                cudaMemcpy(shape , srcshape , sizeof(int64_t) * N, cpyMode);
                cudaMemcpy(stride, srcstride, sizeof(int64_t) * N, cpyMode);
            }
            break;
#endif 
        default:
            Free();
            throw std::runtime_error("Unsupported device when copying (" + std::to_string(device) + "/" + std::to_string(srcDevice) + ")!");

        }
    }
};