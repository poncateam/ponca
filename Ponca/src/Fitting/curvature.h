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


/*!
    \brief Extension to compute curvature values from the Weingarten map \f$ \frac{d N}{d \mathbf{x}} \f$
    \inherit Concept::FittingExtensionConcept

    This class extracts curvature informations from the spatial derivatives of the normal field \f$ N \f$.
    It first assemble a 2x2 matrix representation of the shape operator, and then performs an eigenvalue decomposition
    using Eigen::SelfAdjointEigenSolver::computeDirect.

    The previous basket elements must provide a \c dNormal() method returning a 3x3 matrix.
    If more than one basket element provide a \c dNormal() member, then the last one will be used.

    \warning This class is valid only in 3D.
    \todo Add a compile time check for the working dimension

    \todo Check if PROVIDES_PRINCIPALE_CURVATURES is same as BaseCurvatureEstimator::PROVIDES_PRINCIPALE_CURVATURES

    \ingroup fitting
*/
template < class DataPoint, class _WFunctor, typename T>
class CurvatureEstimator : public T
{
private:
    typedef T Base;

protected:
    enum
    {
        Check = Base::PROVIDES_NORMAL_SPACE_DERIVATIVE,
        PROVIDES_PRINCIPALE_CURVATURES
    };

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using MatrixType = typename Base::MatrixType; /*!< \brief Inherited matrix type*/

private:
    typedef Eigen::Matrix<Scalar,3,2> Mat32; /*!< \brief Matrix type for tangent plane basis */
    typedef Eigen::Matrix<Scalar,2,2> Mat22; /*!< \brief Matrix type for shape operator */

private:
    Scalar m_k1 {0}, m_k2 {0};
    VectorType m_v1 {VectorType::Zero()}, m_v2{VectorType::Zero()};

    static_assert ( DataPoint::Dim == 3, "BaseCurvatureEstimator is only valid in 3D");

public:
    /*! \brief Default constructor */
    PONCA_MULTIARCH inline CurvatureEstimator() = default;

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize();

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
    PONCA_MULTIARCH inline Scalar kMean() const { return (m_k1 + m_k2)/2.;}

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_k1 * m_k2;}

    //! \brief Compute principal curvature directions
    //!
    //! The tangent plane can be calculated from the normal vector or from its
    //! derivatives, depending of the useNormal parameter
    //!
    //! The finalize() method calls this function with useNormal=false by
    //! default.
    //!
    PONCA_MULTIARCH inline FIT_RESULT computeCurvature(bool useNormal = false);

protected:
    //! \brief Compute a tangent plane basis
    /// \todo Uniformize with tangentplane basis
    PONCA_MULTIARCH inline Mat32 tangentPlane(bool useNormal = false) const;
};

#include "curvature.hpp"

} //namespace Ponca
