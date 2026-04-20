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

        /*!
         * \brief Get access to the stored points of evaluation
         * \return Position of the local basis center
         */
        PONCA_MULTIARCH [[nodiscard]] inline const VectorType& evalPos() const { return m_p; }

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
}; // namespace Ponca
