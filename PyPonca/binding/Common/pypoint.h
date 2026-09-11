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

    static constexpr unsigned int Dim = _Dim;

    PONCA_MULTIARCH inline PointNormalBinding(const Scalar* _pos, const Scalar* _normal)
        : m_pos(Eigen::Map<const VectorType>(_pos)), m_normal(Eigen::Map<const VectorType>(_normal))
    {
    }

    PONCA_MULTIARCH [[nodiscard]] inline const Eigen::Map<const VectorType>& pos() const { return m_pos; }
    PONCA_MULTIARCH [[nodiscard]] inline const Eigen::Map<const VectorType>& normal() const { return m_normal; }

private:
    const Eigen::Map<const VectorType> m_pos, m_normal;
};

/**
 * \brief Utility class to provide Ponca-compatible iterators over ndarray
 */
template <typename _Scalar, unsigned int _Dim>
struct PyPointCloud
{
    friend struct iterator;

    using Point      = PointNormalBinding<_Scalar, _Dim>;
    using Scalar     = typename Point::Scalar;
    using VectorType = typename Point::VectorType;

    // Necessary for compatibility with KDTree
    using value_type = Point;

    inline static const std::string PointName = ManglePoint<Point>() + MangledPointName;
    inline static const std::string ClassName = "PointCloud" + PointName;
    static constexpr const std::array<_Scalar, _Dim> NoData{};

    PyPointCloud()
    {
        // We emulate data with an array of stride 0
        void* data = reinterpret_cast<void*>(const_cast<Scalar*>(&NoData[0]));
        m_pos      = PyVectorArray<Point>(data, {1, 1}, nb::handle(), {0, sizeof(Scalar)});
        m_normals  = PyVectorArray<Point>(data, {1, 1}, nb::handle(), {0, sizeof(Scalar)});
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

        if (_pos.stride(1) != 1)
            throw std::runtime_error("PointCloud does not support non contiguous coordinates on second dimensions.");

        // We emulate normal data with an array of stride 0
        void* data = reinterpret_cast<void*>(const_cast<Scalar*>(&NoData[0]));
        m_normals  = PyVectorArray<Point>(data, {_pos.shape(0), 1}, nb::handle(), {0, sizeof(Scalar)});
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

        if (m_pos.stride(1) != 1 || m_normals.stride(1) != 1)
            throw std::runtime_error("PointCloud does not support non contiguous coordinates on second dimensions.");

        if (m_pos.dtype() != m_normals.dtype())
            throw std::runtime_error("Type mismatch between position and normals.");
    }

    auto device_type() const { return m_pos.device_type(); }

    /**
     * \brief Returns the number of points within the array
     */
    size_t size() const { return m_pos.shape(0); }

    /**
     * \brief Accessor
     *
     * No bound checking is performed.
     *
     * \param i Point index
     */
    Point operator[](size_t i) const
    {
        const Scalar* pos    = m_pos.data() + i * m_pos.stride(0);
        const Scalar* normal = m_normals.data() + i * m_normals.stride(0);

        return Point(pos, normal);
    }

    /**
     * \brief Minimal iterator for compatibility with Ponca
     */
    struct iterator
    {
        iterator(const PyPointCloud* _cloud, size_t _idx) : m_cloud(_cloud), m_index(_idx) {}

        // We assume that checking for the cloud adress is useless

        /**
         * \brief Equality comparison
         *
         * \param other The other iterator to check equality with
         */
        bool operator==(const iterator& other) { return m_index == other.m_index; }

        /**
         * \brief Difference comparison
         *
         * \param other The other iterator to check difference with
         */
        bool operator!=(const iterator& other) { return m_index != other.m_index; }

        /**
         * \brief Prefix increment
         */
        iterator& operator++()
        {
            m_index++;
            return *this;
        }

        /**
         * \brief Postfix increment
         */
        iterator operator++(int)
        {
            iterator tmp = *this;
            ++m_index;
            return *this;
        }

        /**
         * \brief Prefix increment
         */
        Point operator*() const { return m_cloud->operator[](m_index); }

    private:
        // TODO: Evaluate if storing two pointers + two strides is faster !
        const PyPointCloud* m_cloud;
        size_t m_index;
    };

    /**
     * \brief Begining iterator
     */
    iterator begin() const { return iterator{this, 0}; }

    /**
     * \brief Begining iterator
     */
    iterator cbegin() const { return iterator{this, 0}; }

    /**
     * \brief End iterator
     */
    iterator end() const { return iterator{this, m_pos.shape(0)}; }

    /**
     * \brief Begining iterator
     */
    iterator cend() const { return iterator{this, m_pos.shape(0)}; }

private:
    PyVectorArray<Point> m_pos;
    PyVectorArray<Point> m_normals;
};

