/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca
{
    /*!
     * \brief NeighborhoodFrame that express 3d points relatively to a prescribed center.
     *
     * This class is useful to get all coordinates centered around a point, which ease weights computation, and
     * limits issues with big numbers and rounding errors.
     *
     * Express points \f$\mathbf{x}\f$ relatively to a center \f$\mathbf{p}\f$, ie.
     * \f$\mathbf{x}'=\mathbf{x}-\mathbf{p}\f$.
     * This frame does not apply rotation.
     *
     * \note This class is designed to serve as a Base class for the definition of NeighborFilters, and cannot be
     * used directly
     *
     * \tparam DataPoint Point type used for computation
     */
    template <class DataPoint>
    class CenteredNeighborhoodFrame
    {
    public:
        /*! \brief Scalar type from DataPoint */
        using Scalar = typename DataPoint::Scalar;
        /*! \brief Vector type from DataPoint */
        using VectorType = typename DataPoint::VectorType;

        /// \brief Flag indicating that this class modifies the coordinates when passing from global to local
        static constexpr bool hasLocalFrame = true;

        PONCA_MULTIARCH inline explicit CenteredNeighborhoodFrame(const VectorType& _evalPos = VectorType::Zero())
            : m_p(_evalPos)
        {
        }

        PONCA_MULTIARCH inline explicit CenteredNeighborhoodFrame(const DataPoint& _evalPoint) : m_p(_evalPoint.pos())
        {
        }

        PONCA_MULTIARCH inline VectorType center() { return m_p; }
        PONCA_MULTIARCH inline VectorType center() const { return m_p; }

        /// \brief Change neighborhood frame (move basis center)
        /// \warning Calling this method invalidates any primitive defined relatively to the previous frame.
        ///          In this situation, it is recommended to use ProvidesImplicitPrimitive::changeBasis instead.
        PONCA_MULTIARCH inline void changeNeighborhoodFrame(const VectorType& _newEvalPos) { m_p = _newEvalPos; };

        /*!
         * \brief Convert query from local to global coordinate system, such as
         * \f$\mathbf{x}=\mathbf{x}'+\mathbf{p}\f$.
         * \param _q Vector expressed relatively to the basis center
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return Vector expressed in the global coordinate system
         *
         * \see convertToLocalBasis
         */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType convertToGlobalBasis(const VectorType& _q,
                                                                             bool _isPositionVector = true) const
        {
            return (_isPositionVector ? (_q + m_p) : _q);
        }

        /*!
         * \brief Convert query from global to local coordinate system, such as
         * \f$\mathbf{x}'=\mathbf{x}-\mathbf{p}\f$.
         *
         * \param _q Input Vector in global coordinate system
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return Vector expressed relatively to the basis center
         *
         * \see convertToGlobalBasis
         */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType convertToLocalBasis(const VectorType& _q,
                                                                            bool _isPositionVector = true) const
        {
            return (_isPositionVector ? (_q - m_p) : _q);
        }

    private:
        VectorType m_p; /*!< \brief basis center */
    };
    /*!
     * \brief NeighborhoodFrame that keep points in the global frame without applying any transformation
     * This class is useful to compute direct fits in the embedding space, without paying the cost to express
     * neighbors relatively to an evaluation point.
     *
     * \warning In case the data have strong magnitude (e.g., georeferenced point clouds), it is recommended to use
     * CenteredNeighborhoodFrame instead.
     *
     * \note This class is designed to serve as a Base class for the definition of NeighborFilters, and cannot be
     * used directly
     *
     * \tparam DataPoint Point type used for computation
     */
    template <class DataPoint>
    class GlobalNeighborhoodFrame
    {
    public:
        /*! \brief Scalar type from DataPoint */
        using Scalar = typename DataPoint::Scalar;
        /*! \brief Vector type from DataPoint */
        using VectorType = typename DataPoint::VectorType;

        /// \brief Flag indicating that this class does not modify the coordinates when passing from global to local
        static constexpr bool hasLocalFrame = false;

        PONCA_MULTIARCH inline explicit GlobalNeighborhoodFrame(const VectorType& /*_evalPos*/ = VectorType::Zero()) {}

        /// \brief Change neighborhood frame (has no effect for global basis)
        /// \warning Calling this method invalidates any primitive defined relatively to the previous frame.
        ///          In this situation, it is recommended to use ProvidesImplicitPrimitive::changeBasis instead.
        PONCA_MULTIARCH inline void changeNeighborhoodFrame(const VectorType& /*_newEvalPos*/) {};

        PONCA_MULTIARCH inline VectorType center() { return VectorType::Zero(); }
        PONCA_MULTIARCH inline VectorType center() const { return VectorType::Zero(); }
        /*!
         * \brief Convert position from local to global coordinate system : does nothing as this is global frame
         * \param _q Position in local coordinate
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return _q
         */
        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& convertToGlobalBasis(
            const VectorType& _q, bool /*_isPositionVector*/ = true) const
        {
            return _q;
        }

        /*!
         * \brief Convert query from global to local coordinate system : does nothing as this is global frame
         * \param _q Query in global coordinate
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return _q
         */
        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& convertToLocalBasis(
            const VectorType& _q, bool /*_isPositionVector*/ = true) const
        {
            return _q;
        }
    };

    /*!
     * \brief NeighborhoodFrame that express points relatively to an ellipsoidal basis.
     *
     * The frame is defined by a center \f$\mathbf{p}\f$ and `DataPoint::Dim` axis vectors
     * \f$\mathbf{a}_0, \dots, \mathbf{a}_{d-1}\f$, stored as the columns of a matrix \f$\mathbf{A}\f$.
     * The direction of each axis gives the orientation of the ellipsoid, and its norm gives the length of the
     * corresponding semi-axis. Points are expressed as \f$\mathbf{x}'=\mathbf{A}^{-1}(\mathbf{x}-\mathbf{p})\f$, so
     * that the unit ball in local coordinates corresponds to the ellipsoid in global coordinates.
     *
     * - With identity axes, this frame is equivalent to CenteredNeighborhoodFrame.
     * - Vectors that are not positions are only transformed by the linear part of the frame
     * (\f$\mathbf{x}'=\mathbf{A}^{-1}\mathbf{x}\f$) (which is not correct for normal vectors...)
     *
     * \note This class is designed to serve as a Base class for the definition of NeighborFilters, and cannot be
     * used directly
     *
     * \tparam DataPoint Point type used for computation
     */
    template <class DataPoint>
    class EllipsoidalNeighborhoodFrame
    {
    public:
        /*! \brief Scalar type from DataPoint */
        using Scalar = typename DataPoint::Scalar;
        /*! \brief Vector type from DataPoint */
        using VectorType = typename DataPoint::VectorType;
        /*! \brief Matrix type from DataPoint */
        using MatrixType = typename DataPoint::MatrixType;

        /// \brief Flag indicating that this class modifies the coordinates when passing from global to local
        static constexpr bool hasLocalFrame = true;

        /// \brief Build the frame from its center and a matrix whose columns are the ellipsoid axes
        PONCA_MULTIARCH inline explicit EllipsoidalNeighborhoodFrame(const VectorType& _evalPos = VectorType::Zero(),
                                                                     const MatrixType& _axes = MatrixType::Identity())
            : m_p(_evalPos)
        {
            setAxes(_axes);
        }

        /// \brief Build the frame from its center and a pointer to `DataPoint::Dim` axis vectors
        PONCA_MULTIARCH inline EllipsoidalNeighborhoodFrame(const VectorType& _evalPos, const VectorType* _axes)
            : m_p(_evalPos)
        {
            setAxes(_axes);
        }

        PONCA_MULTIARCH inline explicit EllipsoidalNeighborhoodFrame(const DataPoint& _evalPoint,
                                                                     const MatrixType& _axes = MatrixType::Identity())
            : EllipsoidalNeighborhoodFrame(_evalPoint.pos(), _axes)
        {
        }

        /// \copydoc EllipsoidalNeighborhoodFrame(const VectorType&, const VectorType*)
        PONCA_MULTIARCH inline EllipsoidalNeighborhoodFrame(const DataPoint& _evalPoint, const VectorType* _axes)
            : EllipsoidalNeighborhoodFrame(_evalPoint.pos(), _axes)
        {
        }

        PONCA_MULTIARCH inline VectorType center() { return m_p; }
        PONCA_MULTIARCH inline VectorType center() const { return m_p; }

        /// \brief Matrix whose columns are the ellipsoid axes
        PONCA_MULTIARCH inline const MatrixType& axes() const { return m_axes; }
        /// \brief Inverse of the matrix whose columns are the ellipsoid axes
        PONCA_MULTIARCH inline const MatrixType& inverseAxes() const { return m_invAxes; }
        /// \brief i-th axis of the ellipsoid
        PONCA_MULTIARCH inline VectorType axis(int _i) const { return m_axes.col(_i); }

        /// \brief Set the ellipsoid axes from a matrix whose columns are the axes
        PONCA_MULTIARCH inline void setAxes(const MatrixType& _axes)
        {
            m_axes    = _axes;
            m_invAxes = _axes.inverse();
        }

        /// \brief Set the ellipsoid axes from a pointer to `DataPoint::Dim` contiguous vectors
        PONCA_MULTIARCH inline void setAxes(const VectorType* _axes)
        {
            MatrixType axes;
            for (int i = 0; i < DataPoint::Dim; ++i)
                axes.col(i) = _axes[i];
            setAxes(axes);
        }

        /// \brief Change neighborhood frame (move basis center, axes are left unchanged)
        PONCA_MULTIARCH inline void changeNeighborhoodFrame(const VectorType& _newEvalPos) { m_p = _newEvalPos; };

        /*!
         * \brief Convert query from local to global coordinate system, such as
         * \f$\mathbf{x}=\mathbf{A}\mathbf{x}'+\mathbf{p}\f$.
         * \param _q Vector expressed in the ellipsoidal basis
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return Vector expressed in the global coordinate system
         *
         * \see convertToLocalBasis
         */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType convertToGlobalBasis(const VectorType& _q,
                                                                             bool _isPositionVector = true) const
        {
            return _isPositionVector ? VectorType(m_axes * _q + m_p) : VectorType(m_axes * _q);
        }

        /*!
         * \brief Convert query from global to local coordinate system, such as
         * \f$\mathbf{x}'=\mathbf{A}^{-1}(\mathbf{x}-\mathbf{p})\f$.
         *
         * \param _q Input Vector in global coordinate system
         * \param _isPositionVector Indicate if the input vector `_q` is a position that is influenced by
         * translations (e.g., in contrast to displacement or normal vectors)
         * \return Vector expressed in the ellipsoidal basis
         *
         * \see convertToGlobalBasis
         */
        PONCA_MULTIARCH [[nodiscard]] inline VectorType convertToLocalBasis(const VectorType& _q,
                                                                            bool _isPositionVector = true) const
        {
            return _isPositionVector ? VectorType(m_invAxes * (_q - m_p)) : VectorType(m_invAxes * _q);
        }

    private:
        VectorType m_p;       /*!< \brief basis center */
        MatrixType m_axes;    /*!< \brief ellipsoid axes, stored as columns */
        MatrixType m_invAxes; /*!< \brief inverse of m_axes, cached */
    };
}; // namespace Ponca
