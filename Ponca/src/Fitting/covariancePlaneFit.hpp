/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
CovariancePlaneFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    if (Base::finalize() == STABLE) {
        if (Base::plane().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
        Base::setPlane(Base::m_solver.eigenvectors().col(0), Base::barycenterLocal());
    }

    return Base::m_eCurrentState;
}

template < class DataPoint, class _NFilter, typename T>
typename CovariancePlaneFitImpl<DataPoint, _NFilter, T>::VectorType
CovariancePlaneFitImpl<DataPoint, _NFilter, T>::worldToTangentPlane (const VectorType& _q, bool _isPositionVector) const
{
    return Base::m_solver.eigenvectors().transpose() * Base::getNeighborFilter().convertToLocalBasis(_q,
                                                                                                     _isPositionVector);
}

template < class DataPoint, class _NFilter, typename T>
typename CovariancePlaneFitImpl<DataPoint, _NFilter, T>::VectorType
CovariancePlaneFitImpl<DataPoint, _NFilter, T>::tangentPlaneToWorld (const VectorType& _lq, bool _isPositionVector) const
{
    return Base::getNeighborFilter().convertToGlobalBasis(Base::m_solver.eigenvectors().transpose().inverse() * _lq,
                                                          _isPositionVector);
}



template < class DataPoint, class _NFilter, int DiffType, typename T>
FIT_RESULT
CovariancePlaneDerImpl<DataPoint, _NFilter, DiffType, T>::finalize()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);
    PONCA_MULTIARCH_STD_MATH(numeric_limits);

    Base::finalize();
    // Test if base finalize end on a viable case (stable / unstable)
    if (this->isReady())
    {
      VectorType   barycenter = Base::barycenterLocal();
      VectorArray dBarycenter = Base::barycenterDerivatives();

      // pre-compute shifted eigenvalues to apply the pseudo inverse of C - lambda_0 I
      Scalar epsilon          = Scalar(2) * Eigen::NumTraits<Scalar>::epsilon();
      Scalar consider_as_zero = Scalar(2) * numeric_limits<Scalar>::denorm_min();

      // This is where the limitation to 3d comes from.
      // \fixme Replace shift in 2d subspace by any subspace with co-dimension 1
      Eigen::Matrix<Scalar,2,1> shifted_eivals = Base::m_solver.eigenvalues().template tail<2>().array() - Base::m_solver.eigenvalues()(0);
      if(shifted_eivals(0) < consider_as_zero || shifted_eivals(0) < epsilon * shifted_eivals(1)) shifted_eivals(0) = 0;
      if(shifted_eivals(1) < consider_as_zero) shifted_eivals(1) = 0;


      for(int k=0; k<Base::NbDerivatives; ++k)
      {
        VectorType normal = Base::primitiveGradient();
        // The derivative of 'normal' is the derivative of the smallest eigenvector.
        // Since the covariance matrix is real and symmetric, it is equal to:
        //    n' = - (C - lambda_0 I)^+ C' n
        // Where ^+ denotes the pseudo-inverse.
        // Since we already performed the eigenvalue decomposition of the matrix C,
        // we can directly apply the pseudo inverse by observing that:
        //    (C - lambda_0 I) = V (L - lambda_0 I) V^T
        // where V is the eigenvector matrix, and L the eigenvalue diagonal matrix.
        Eigen::Matrix<Scalar,2,1> z = - Base::m_solver.eigenvectors().template rightCols<2>().transpose() * (Base::m_dCov[k] * normal);
        if(shifted_eivals(0)>0) z(0) /= shifted_eivals(0);
        if(shifted_eivals(1)>0) z(1) /= shifted_eivals(1);
        m_dNormal.col(k) = Base::m_solver.eigenvectors().template rightCols<2>() * z;

        VectorType dDiff = dBarycenter.col(k);
        if(k>0 || !Base::isScaleDer())
          dDiff(Base::isScaleDer() ? k-1 : k) += 1;
        m_dDist(k) = m_dNormal.col(k).dot(barycenter) + normal.dot(dDiff);

        // \fixme we shouldn't need this normalization, however currently the derivatives are overestimated by a factor 2
        m_dNormal /= Scalar(2.);
      }
    }

    return Base::m_eCurrentState;
}
