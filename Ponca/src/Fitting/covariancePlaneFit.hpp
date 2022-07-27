/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
CovariancePlaneFitImpl<DataPoint, _WFunctor, T>::finalize ()
{
    if (Base::finalize() == STABLE)
        Base::setPlane(Base::m_solver.eigenvectors().col(0), Base::barycenter());

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename CovariancePlaneFitImpl<DataPoint, _WFunctor, T>::VectorType
CovariancePlaneFitImpl<DataPoint, _WFunctor, T>::worldToTangentPlane (const VectorType& _q) const
{
  if (ignoreTranslation)
    return Base::m_solver.eigenvectors().transpose() * _q;
  else {
    // apply rotation and translation to get uv coordinates
    return Base::m_solver.eigenvectors().transpose() * (Base::m_w.convertToLocalBasis(_q));
  }
}

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename CovariancePlaneFitImpl<DataPoint, _WFunctor, T>::VectorType
CovariancePlaneFitImpl<DataPoint, _WFunctor, T>::tangentPlaneToWorld (const VectorType& _lq) const
{
  if (ignoreTranslation)
    return Base::m_solver.eigenvectors().transpose().inverse() * _lq;
  else {
    return Base::m_solver.eigenvectors().transpose().inverse() * _lq + Base::m_w.basisCenter();
  }
}



template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
CovariancePlaneDerImpl<DataPoint, _WFunctor, DiffType, T>::finalize()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    Base::finalize();
    // Test if base finalize end on a viable case (stable / unstable)
    if (this->isReady())
    {
      VectorType   barycenter = Base::barycenter();
      VectorArray dBarycenter = Base::barycenterDerivatives();

      // pre-compute shifted eigenvalues to apply the pseudo inverse of C - lambda_0 I
      Scalar epsilon          = Scalar(2) * Eigen::NumTraits<Scalar>::epsilon();
      Scalar consider_as_zero = Scalar(2) * std::numeric_limits<Scalar>::denorm_min();

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
