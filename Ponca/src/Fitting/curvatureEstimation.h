/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

namespace Ponca
{

template < class DataPoint, class _WFunctor, typename T>
/**
 *
 * \brief Base class for any 3d curvature estimator: holds k1, k2 and associated vectors
 *
 * \todo Check if PROVIDES_PRINCIPAL_CURVATURES is same as CurvatureEstimator::PROVIDES_PRINCIPAL_CURVATURES
 *
 * \ingroup fitting
 */
class BaseCurvatureEstimator : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { PROVIDES_PRINCIPAL_CURVATURES };

protected:
    /// \brief Principal curvature with highest absolute magnitude
    Scalar m_k1 {0},
    /// \brief Principal curvature with smallest absolute magnitude
           m_k2 {0};
    /// \brief Direction associated to the principal curvature with highest absolute magnitude
    VectorType m_v1 {VectorType::Zero()},
    /// \brief Direction associated to the principal curvature with highest smallest magnitude
               m_v2 {VectorType::Zero()};

    static_assert ( DataPoint::Dim == 3, "BaseCurvatureEstimator is only valid in 3D");

public:
    PONCA_EXPLICIT_CAST_OPERATORS(BaseCurvatureEstimator,baseCurvatureEstimator)
    PONCA_FITTING_DECLARE_INIT

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
};


/*!
 * \brief Extension to compute curvature values based on a covariance analysis
 * of normal vectors of neighbors.
 *
 * A 3D covariance matrix is computed from the normals of the neighbors and the
 * two principal curvature values and directions are given by the two extreme
 * eigenvalues and associated eigenvectors of the covariance matrix
 * \cite Liang:1990:RRSS.
 *
 * \todo Remove direct inheritance (use PROVIDE system instead)
 *
 * \warning This class is valid only in 3D.
 * \ingroup fitting
 */
template < class DataPoint, class _WFunctor, typename T>
class NormalCovarianceCurvature : public BaseCurvatureEstimator<DataPoint,_WFunctor,T>
{
private:
    // \todo Remove this, replace by PONCA_FITTING_DECLARE_DEFAULT_TYPES
    typedef BaseCurvatureEstimator<DataPoint,_WFunctor,T> Base;

    //TODO(thib) check the curvature values that might be wrong
    //TODO(thib) use weighting function
    //TODO(thib) which eigenvectors should be selected ? extreme of maximal ?

public:
    typedef typename Base::Scalar          Scalar;      /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType      VectorType;  /*!< \brief Inherited vector type*/
    typedef typename DataPoint::MatrixType MatrixType;  /*!< \brief Matrix type inherited from DataPoint*/
    /*! \brief Solver used to analyse the covariance matrix*/
    typedef Eigen::SelfAdjointEigenSolver<MatrixType> Solver;

protected:
    MatrixType m_cov;   /*!< \brief Covariance matrix of the normal vectors */
    VectorType m_cog;   /*!< \brief Gravity center of the normal vectors \fixme Use MeanNormal */
    Solver m_solver;    /*!< \brief Solver used to analyse the covariance matrix */

public:
    PONCA_EXPLICIT_CAST_OPERATORS(NormalCovarianceCurvature,normalCovarianceCurvature)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
};


/*!
 * \brief Extension to compute curvature values based on a covariance analysis
 * of normal vectors of neighbors projected onto the tangent plane.
 *
 * A 2D covariance matrix is computed from the projections of normals of the
 * neighbors and the two principal curvature values and directions are given by
 * the eigenvalues and associated eigenvectors of the covariance matrix
 * \cite Berkmann:1994:CSG.
 *
 * \note This procedure requires two passes, the first one for plane fitting
 * and local frame estimation, and the second one for covariance analysis.
 * \warning This class is valid only in 3D.
 * \ingroup fitting
 */
template < class DataPoint, class _WFunctor, typename T>
class ProjectedNormalCovarianceCurvature : public BaseCurvatureEstimator<DataPoint,_WFunctor,T>
{
private:
    // \todo Remove this, replace by PONCA_FITTING_DECLARE_DEFAULT_TYPES
    typedef BaseCurvatureEstimator<DataPoint,_WFunctor,T> Base;

    //TODO(thib) check the curvature values that might be wrong
    //TODO(thib) use weighting function

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE
    };

    /// \fixme Use same pass management than MongePatch
    enum PASS : int
    {
        FIRST_PASS = 0,
        SECOND_PASS,
        PASS_COUNT
    };

public:
    typedef typename Base::Scalar          Scalar;      /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType      VectorType;  /*!< \brief Inherited vector type*/
    typedef typename DataPoint::MatrixType MatrixType;  /*!< \brief Matrix type inherited from DataPoint*/
    //TODO(thib) use of Eigen::RowAtCompileTime-1 ?
    typedef Eigen::Matrix<Scalar,2,2> Mat22;
    typedef Eigen::Matrix<Scalar,3,2> Mat32;
    typedef Eigen::Matrix<Scalar,2,1> Vector2;
    typedef typename VectorType::Index Index;
    /*! \brief Solver used to analyse the covariance matrix*/
    typedef Eigen::SelfAdjointEigenSolver<Mat22> Solver;

protected:
    Vector2 m_cog;      /*!< \brief Gravity center */
    Mat22 m_cov;        /*!< \brief Covariance matrix */
    Solver m_solver;    /*!< \brief Solver used to analyse the covariance matrix */
    PASS m_pass;        /*!< \brief Current pass */
    Mat32 m_tframe;     /*!< \brief Tangent frame */

public:
    PONCA_EXPLICIT_CAST_OPERATORS(ProjectedNormalCovarianceCurvature,projectedNormalCovarianceCurvature)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
};

#include "curvatureEstimation.hpp"

} //namespace Ponca
