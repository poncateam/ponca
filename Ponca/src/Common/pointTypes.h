/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

/*!
  \file pointTypes.h
  \brief Useful user-end DataPoint types
*/

namespace Ponca {
    //! [PointPositionNormal]
    /*! \brief Point with position and normal vector
     */
    template<typename _Scalar, int _Dim>
    class PointPositionNormal
    {
    public:
        enum {Dim = _Dim};
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Dim,   1>	VectorType;
        typedef Eigen::Matrix<Scalar, Dim, Dim>	MatrixType;

        PONCA_MULTIARCH inline PointPositionNormal(
                const VectorType &pos    = VectorType::Zero(),
                const VectorType& normal = VectorType::Zero()
        ) : m_pos(pos), m_normal(normal) {}

        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& pos()    const { return m_pos; }
        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& normal() const { return m_normal; }
        PONCA_MULTIARCH [[nodiscard]] inline VectorType& pos()    { return m_pos; }
        PONCA_MULTIARCH [[nodiscard]] inline VectorType& normal() { return m_normal; }

    private:
        VectorType m_pos, m_normal;
    };
    //! [PointPositionNormal]

    //! [PointPosition]
    /// Point with position, without attribute
    template<typename _Scalar, int _Dim>
    class PointPosition
    {
    public:
        enum {Dim = _Dim};
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Dim,   1>	VectorType;
        typedef Eigen::Matrix<Scalar, Dim, Dim>	MatrixType;

        PONCA_MULTIARCH inline PointPosition(
            const VectorType &pos = VectorType::Zero()
        ) : m_pos(pos) {}

        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& pos()    const { return m_pos; }
        PONCA_MULTIARCH [[nodiscard]] inline VectorType& pos()    { return m_pos; }

    private:
        VectorType m_pos;
    };
    //! [PointPosition]


    /*! \brief Variant of the \ref PointPositionNormal data type that uses external raw data.
     * Using this approach, ones can use the ponca library with already existing
     * data-structures and without any data-duplication.
     *
     * In this example, we use this class to map an interlaced raw array containing
     * both point normals and coordinates.
     *
     * \see PointPositionNormal
     */
    template<typename _Scalar, int _Dim>
    class PointPositionNormalBinding
    {
    public:
        enum {Dim = _Dim};
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
        typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

        PONCA_MULTIARCH inline PointPositionNormalBinding(
            const Scalar* _interlacedArray, int _pId
        ) : m_pos   (Eigen::Map< const VectorType >(_interlacedArray + Dim*2*_pId  )),
            m_normal(Eigen::Map< const VectorType >(_interlacedArray + Dim*2*_pId+Dim))
        {}

        PONCA_MULTIARCH [[nodiscard]] inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
        PONCA_MULTIARCH [[nodiscard]] inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

    private:
        Eigen::Map< const VectorType > m_pos, m_normal;
    };
}
