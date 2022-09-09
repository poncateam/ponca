/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Core>

namespace Ponca
{
    template < class DataPoint, class _WFunctor, int DiffType, typename T>
/**
 *
 * \brief Base class for any 3d curvature estimator: holds \f$k_1\f$, \f$k_2\f$ and associated vectors
 */
    class CurvatureEstimatorBase : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

    protected:
        enum
        {
            PROVIDES_PRINCIPAL_CURVATURES
        };

    private:
        /// \brief Principal curvature with highest absolute magnitude
        Scalar m_k1 {0},
        /// \brief Principal curvature with smallest absolute magnitude
        m_k2 {0};
        /// \brief Direction associated to the principal curvature with highest absolute magnitude
        VectorType m_v1 {VectorType::Zero()},
        /// \brief Direction associated to the principal curvature with highest smallest magnitude
        m_v2 {VectorType::Zero()};

        /// \brief Internal state indicating if the curvature are set to default (false), or have been computed (true)
        bool m_isValid {false};

        static_assert ( DataPoint::Dim == 3, "CurvatureEstimatorBase is only valid in 3D");

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(CurvatureEstimatorBase, curvatureEstimatorBase)
        PONCA_FITTING_DECLARE_INIT

        /// \brief Returns true if contains valid curvature values (and not default ones)
        PONCA_MULTIARCH inline bool isValid() const { return m_isValid; }

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
        //! \brief Returns an estimate of the first principal curvature value
        //!
        //! It is the greatest curvature in <b>absolute value</b>.
        PONCA_MULTIARCH inline Scalar k1() const { return m_k1; }

        //! \brief Returns an estimate of the second principal curvature value
        //!
        //! It is the smallest curvature in <b>absolute value</b>.
        PONCA_MULTIARCH inline Scalar k2() const { return m_k2; }

        //! \brief Returns an estimate of the first principal curvature direction
        //!
        //! It is the greatest curvature in <b>absolute value</b>.
        PONCA_MULTIARCH inline VectorType k1Direction() const { return m_v1; }

        //! \brief Returns an estimate of the second principal curvature direction
        //!
        //! It is the smallest curvature in <b>absolute value</b>.
        PONCA_MULTIARCH inline VectorType k2Direction() const { return m_v2; }

        //! \brief Returns an estimate of the mean curvature
        PONCA_MULTIARCH inline Scalar kMean() const { return (m_k1 + m_k2)/Scalar(2);}

        //! \brief Returns an estimate of the Gaussian curvature
        PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_k1 * m_k2;}

    protected:
        /// \brief Set curvature values. To be called in finalize() by child classes
        ///
        /// Principal curvatures are re-ordered depending on their magnitude,
        /// such that \f$ |k_1| > |k_2| \f$
        PONCA_MULTIARCH inline void setCurvatureValues(Scalar k1, Scalar k2, const VectorType& v1, const VectorType& v2);
    };

} //namespace Ponca

#include "curvature.hpp"
