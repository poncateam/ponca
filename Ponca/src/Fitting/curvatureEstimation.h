/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

namespace Ponca
{
/*!
    \brief Extension to compute curvature values from the Weingarten map \f$ \frac{d N}{d \mathbf{x}} \f$
    \inherit Concept::FittingExtensionConcept

    This class extracts curvature information from the spatial derivatives of the normal field \f$ N \f$.
    It first assemble a 2x2 matrix representation of the shape operator, and then performs an eigenvalue decomposition
    using Eigen::SelfAdjointEigenSolver::computeDirect.
*/
    template < class DataPoint, class _NFilter, int DiffType, typename T>
    class NormalDerivativesCurvatureEstimator : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

    protected:
        enum
        {
            Check = Base::PROVIDES_NORMAL_DERIVATIVE && Base::PROVIDES_PRINCIPAL_CURVATURES
        };

    private:
        typedef Eigen::Matrix<Scalar,3,2> Mat32; /*!< \brief Matrix type for tangent plane basis \fixme formalize tangent plane basis */
        typedef Eigen::Matrix<Scalar,2,2> Mat22; /*!< \brief Matrix type for shape operator */

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(NormalDerivativesCurvatureEstimator,normalDerivativesCurvatureEstimator)
        PONCA_FITTING_DECLARE_FINALIZE

    private:
        //! \brief Compute principal curvature directions relatively to the tangent plane
        //! \see tangentPlane
        //! The finalize() method calls this function with useNormal=false by default.
        //! \todo Add a way to give user control to the tangent plane estimation
        //!
        PONCA_MULTIARCH inline FIT_RESULT computeCurvature(bool useNormal = false);

    protected:
        //! \brief Compute a tangent plane basis
        //!
        //! The tangent plane can be calculated from the normal vector or from its
        //! derivatives, depending of the useNormal parameter
        //! \todo Uniformize with tangentplane basis: these computations are not part of NormalDerivativesCurvature
        PONCA_MULTIARCH inline Mat32 tangentPlane(bool useNormal = false) const;
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
 * \todo Refactor curvature estimators, and link to tangent plane
 *
 * \warning Not it test suite, to be added !
 *
 * \warning This class is valid only in 3D.
 */
template < class DataPoint, class _NFilter, int DiffType, typename T>
class NormalCovarianceCurvatureEstimator : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE
PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_PRINCIPAL_CURVATURES
    };

    //TODO(thib) check the curvature values that might be wrong
    //TODO(thib) use weighting function
    //TODO(thib) which eigenvectors should be selected ? extreme of maximal ?

public:
    /*! \brief Solver used to analyse the covariance matrix*/
    typedef Eigen::SelfAdjointEigenSolver<MatrixType> Solver;

protected:
    MatrixType m_cov;   /*!< \brief Covariance matrix of the normal vectors \todo We have this somewhere else */
    VectorType m_cog;   /*!< \brief Gravity center of the normal vectors \todo Use MeanNormal */
    Solver m_solver;    /*!< \brief Solver used to analyse the covariance matrix */

public:
    PONCA_EXPLICIT_CAST_OPERATORS_DER(NormalCovarianceCurvatureEstimator, normalCovarianceCurvatureEstimator)
    PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE
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
 *
 * \warning Not it test suite, to be added !
 */
template < class DataPoint, class _NFilter, int DiffType, typename T>
class ProjectedNormalCovarianceCurvatureEstimator : public T
{
    //TODO(thib) check the curvature values that might be wrong
    //TODO(thib) use weighting function

PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE
PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_PRINCIPAL_CURVATURES &&
                Base::PROVIDES_PLANE // \todo This class relies on the primitiveGradient, so update this
    };

    /// \todo Use same pass management than MongePatch
    enum PASS : int
    {
        FIRST_PASS = 0,
        SECOND_PASS,
        PASS_COUNT
    };

public:
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
    PONCA_EXPLICIT_CAST_OPERATORS_DER(ProjectedNormalCovarianceCurvatureEstimator, projectedNormalCovarianceCurvature)
    PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE
};

#include "curvatureEstimation.hpp"

} //namespace Ponca
