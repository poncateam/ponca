/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <type_traits>
#include <string>

#include "DeviceUtils.h"
#include "Mangling.h"
#include "pytypes.h"

namespace nb = nanobind;

inline static constexpr const char* MangledPointName = "PN";

/**
 * \brief Creates a Ponca point from two pointers
 */
template <typename _Scalar, unsigned int _Dim>
struct PointNormalBinding
{
    using Scalar     = _Scalar;
    using VectorType = Eigen::Matrix<Scalar, _Dim, 1>;
    using MatrixType = Eigen::Matrix<Scalar, _Dim, _Dim>;

    using Stride = Eigen::InnerStride<Eigen::Dynamic>;
    using Map    = Eigen::Map<const VectorType, Eigen::Unaligned, Stride>;

    static constexpr unsigned int Dim = _Dim;

    PONCA_MULTIARCH inline PointNormalBinding(const Scalar* _pos, size_t _posStride, const Scalar* _normal,
                                              size_t _normalStride)
        : m_pos(_pos, Stride(_posStride)), m_normal(_normal, Stride(_normalStride))
    {
    }

    PONCA_MULTIARCH [[nodiscard]] inline const Map& pos() const { return m_pos; }
    PONCA_MULTIARCH [[nodiscard]] inline const Map& normal() const { return m_normal; }

private:
    const Map m_pos, m_normal;
};

/**
 * \brief Ponca-compatible container of points over raw arrays
 *
 * This is the representation of a point cloud used by computations, on every device:
 * it only holds raw pointers and strides (and is trivially copyable), so it can be used
 * in device code, given the arrays are stored on the device.
 *
 * It does not own any data, see PyPointCloud.
 */
template <typename _Scalar, unsigned int _Dim>
struct PyPointCloudView
{
    using Point  = PointNormalBinding<_Scalar, _Dim>;
    using Scalar = typename Point::Scalar;

    // Necessary for compatibility with KDTree
    using value_type = Point;

    /**
     * \brief Ctor
     *
     * Strides are given in elements, for each dimension of the (N, Dim) arrays.
     *
     * \param pos Pointer to the first position
     * \param posStride Strides of the position array
     * \param normal Pointer to the first normal
     * \param normalStride Strides of the normal array
     * \param size Number of points
     */
    PONCA_MULTIARCH PyPointCloudView(const Scalar* pos, const int64_t posStride[2], const Scalar* normal,
                                     const int64_t normalStride[2], size_t size)
        : m_pos(pos), m_normal(normal), m_posStride{posStride[0], posStride[1]},
          m_normalStride{normalStride[0], normalStride[1]}, m_size(size)
    {
    }

    /**
     * \brief Returns the number of points
     */
    PONCA_MULTIARCH size_t size() const { return m_size; }

    /**
     * \brief Accessor
     *
     * No bound checking is performed.
     *
     * \param i Point index
     */
    PONCA_MULTIARCH Point operator[](size_t i) const
    {
        return Point(m_pos + i * m_posStride[0], m_posStride[1], m_normal + i * m_normalStride[0], m_normalStride[1]);
    }

    /**
     * \brief Iterator for compatibility with Ponca and STL
     */
    struct iterator
    {
        using value_type        = Point;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type;
        using reference         = value_type;
        using iterator_category = std::random_access_iterator_tag;
        using iterator_concept  = std::random_access_iterator_tag;

        PONCA_MULTIARCH iterator(const PyPointCloudView* _cloud, size_t _idx) : m_cloud(_cloud), m_index(_idx) {}

        // We assume that checking for the cloud adress is useless

        /**
         * \brief Equality comparison
         *
         * \param other The other iterator to check equality with
         */
        PONCA_MULTIARCH bool operator==(const iterator& other) const { return m_index == other.m_index; }

        /**
         * \brief Difference comparison
         *
         * \param other The other iterator to check difference with
         */
        PONCA_MULTIARCH bool operator!=(const iterator& other) const { return m_index != other.m_index; }

        /**
         * \brief Prefix increment
         */
        PONCA_MULTIARCH iterator& operator++()
        {
            m_index++;
            return *this;
        }

        /**
         * \brief Prefix decrement
         */
        PONCA_MULTIARCH iterator& operator--()
        {
            m_index--;
            return *this;
        }

        /**
         * \brief Postfix increment
         */
        PONCA_MULTIARCH iterator operator++(int)
        {
            iterator tmp = *this;
            ++m_index;
            return tmp;
        }

        /**
         * \brief Postfix decrement
         */
        PONCA_MULTIARCH iterator operator--(int)
        {
            iterator tmp = *this;
            --m_index;
            return tmp;
        }

        /**
         * \brief General increment
         *
         * \param n Increment amount
         */
        PONCA_MULTIARCH iterator& operator+=(difference_type n)
        {
            m_index += n;
            return *this;
        }

        /**
         * \brief General decrement
         *
         * \param n Decrement amount
         */
        PONCA_MULTIARCH iterator& operator-=(difference_type n)
        {
            m_index -= n;
            return *this;
        }

        PONCA_MULTIARCH difference_type operator-(const iterator& other) const { return m_index - other.m_index; }

        /**
         * \brief Point access
         */
        PONCA_MULTIARCH Point operator*() const { return (*m_cloud)[m_index]; }

        /**
         * \brief Point access
         */
        PONCA_MULTIARCH Point operator->() const { return (*m_cloud)[m_index]; }

    private:
        // TODO: Evaluate if storing two pointers + two strides is faster !
        const PyPointCloudView* m_cloud;
        size_t m_index;
    };

    /**
     * \brief Begining iterator
     */
    PONCA_MULTIARCH iterator begin() const { return iterator{this, 0}; }

    /**
     * \brief Begining iterator
     */
    PONCA_MULTIARCH iterator cbegin() const { return iterator{this, 0}; }

    /**
     * \brief End iterator
     */
    PONCA_MULTIARCH iterator end() const { return iterator{this, m_size}; }

    /**
     * \brief End iterator
     */
    PONCA_MULTIARCH iterator cend() const { return iterator{this, m_size}; }

private:
    const Scalar* m_pos;
    const Scalar* m_normal;
    int64_t m_posStride[2];
    int64_t m_normalStride[2];
    size_t m_size;
};

/**
 * \brief Point cloud given from python
 *
 * This class keeps the arrays alive and checks them. Points are accessed through
 * the view returned by view().
 */
template <typename _Scalar, unsigned int _Dim>
struct PyPointCloud
{
    using Point      = PointNormalBinding<_Scalar, _Dim>;
    using Scalar     = typename Point::Scalar;
    using VectorType = typename Point::VectorType;
    using View       = PyPointCloudView<_Scalar, _Dim>;

    inline static const std::string PointName = ManglePoint<Point>() + MangledPointName;
    inline static const std::string ClassName = "PointCloud" + PointName;
    static constexpr const std::array<_Scalar, _Dim> NoData{};

    PyPointCloud()
    {
        // We emulate data with an array of stride 0 (strides are given in elements)
        void* data = reinterpret_cast<void*>(const_cast<Scalar*>(&NoData[0]));
        m_pos      = PyVectorArray<Point>(data, {1, _Dim}, nb::handle(), {0, 1});
        m_normals  = PyVectorArray<Point>(data, {1, _Dim}, nb::handle(), {0, 1});
    }

    /**
     * \brief Ctor
     *
     * This constructors fills the normal with a 0-stride array
     *
     * \param _pos The position array
     */
    PyPointCloud(const PyVectorArray<Point>& _pos) : m_pos(_pos)
    {
        if (m_pos.ndim() != 2)
            throw std::runtime_error("PointCloud only supports 2D arrays");

        // We emulate normal data with an array of stride 0 (strides are given in elements)
        void* data = reinterpret_cast<void*>(const_cast<Scalar*>(&NoData[0]));
        nb::object owner;
#ifdef __CUDACC__
        // NoData lives on the host, device code needs its own zeros.
        if (_pos.device_type() == nb::device::cuda::value)
        {
            CudaDeviceGuard guard(_pos.device_id());

            Scalar* zeros = nullptr;
            CheckCuda(cudaMalloc(&zeros, _Dim * sizeof(Scalar)), "Device allocation failed");
            owner = nb::capsule(zeros, [](void* p) noexcept { cudaFree(p); });
            CheckCuda(cudaMemset(zeros, 0, _Dim * sizeof(Scalar)), "Device memset failed");

            data = zeros;
        }
#endif
        m_normals = PyVectorArray<Point>(data, {_pos.shape(0), _Dim}, owner, {0, 1}, nb::dtype<Scalar>(),
                                         _pos.device_type(), _pos.device_id());
    }

    /**
     * \brief Ctor
     *
     * \param _pos The position array
     * \param _normals The position array
     */
    PyPointCloud(const PyVectorArray<Point>& _pos, const PyVectorArray<Point>& _normals)
        : m_pos(_pos), m_normals(_normals)
    {
        if (m_pos.ndim() != 2 || m_normals.ndim() != 2)
            throw std::runtime_error("PointCloud only supports 2D arrays");

        if (m_pos.dtype() != m_normals.dtype())
            throw std::runtime_error("Type mismatch between position and normals.");
    }

    auto device_type() const { return m_pos.device_type(); }

    auto device_id() const { return m_pos.device_id(); }

    /**
     * \brief Returns the number of points within the array
     */
    size_t size() const { return m_pos.shape(0); }

    /**
     * \brief Returns the position array
     */
    const PyVectorArray<Point>& positions() const { return m_pos; }

    /**
     * \brief Returns a view over the cloud, usable on the device the arrays are stored on
     */
    View view() const
    {
        const int64_t posStride[2]    = {m_pos.stride(0), m_pos.stride(1)};
        const int64_t normalStride[2] = {m_normals.stride(0), m_normals.stride(1)};
        return View(m_pos.data(), posStride, m_normals.data(), normalStride, size());
    }

private:
    PyVectorArray<Point> m_pos;
    PyVectorArray<Point> m_normals;
};
