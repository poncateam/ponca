/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
 GLSDer <DataPoint, _WFunctor, T>::dtau() const
{
    MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;
    ScalarArray dfield = Base::m_dUc;
    // Recall that tau is the field function at the evaluation point, we thus must take care about
    // its variation when differentiating in space:
    if(this->isScaleDer())
      dfield.template tail<DataPoint::Dim>() += Base::m_ul;

    return (dfield * prattNorm - Base::m_uc * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta() const
{
  // Recall that eta() returns the normal at the evaluation point, therefore,
  // we must take into account the variation of the evaluation point when differentiating wrt space
  // i.e., eta(x) = grad/|grad|, with grad(x) = ul + 2 uq * x, and diff_x(grad) = dul + 2 uq I
  VectorArray dgrad = Base::m_dUl;
  if(this->isSpaceDer())
    dgrad.template rightCols<DataPoint::Dim>().diagonal().array() += Scalar(2)*Base::m_uq;
  Scalar norm  = Base::m_ul.norm();
  Scalar norm3 = norm*norm*norm;
  return dgrad / norm - Base::m_ul * (Base::m_ul.transpose() * dgrad) / norm3;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa() const
{
    MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;

    return Scalar(2.) * (Base::m_dUq * prattNorm - Base::m_uq * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau_normalized() const
{
    return dtau();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta_normalized() const
{
    return Base::m_t * deta();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa_normalized() const
{
    return dkappa() * Base::m_t * Base::m_t;
}




template < class DataPoint, class _WFunctor, typename T>
typename GLSGeomVar <DataPoint, _WFunctor, T>::Scalar
GLSGeomVar <DataPoint, _WFunctor, T>::geomVar(  Scalar wtau, 
                                                Scalar weta,
                                                Scalar wkappa ) const
{
    Scalar dtau   = Base::dtau_normalized().col(0)(0);
    Scalar deta   = Base::deta_normalized().col(0).norm();
    Scalar dkappa = Base::dkappa_normalized().col(0)(0);

    return wtau*dtau*dtau + weta*deta*deta + wkappa*dkappa*dkappa;
}

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
GLSCurvatureHelper<DataPoint, _WFunctor, T>::finalize()
{
    typedef typename VectorType::Index Index;
    typedef Eigen::Matrix<Scalar,3,2> Mat32;
    typedef Eigen::Matrix<Scalar,2,2> Mat22;
    
    MULTIARCH_STD_MATH(sqrt);

    FIT_RESULT bResult = Base::finalize();

    if(bResult != UNDEFINED)
    {
        // Get the object space Weingarten map dN
        MatrixType dN = Base::deta().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);
        
        // Make sure dN is orthogonal to the normal: (optional, does not seem to improve accuracy)
//         VectorType n = Base::eta().normalized();
//         dN = dN - n * n.transpose() * dN;
        
        // Make sure that dN is symmetric:
        // FIXME check why dN is not already symmetric (i.e., round-off errors or error in derivative formulas?)
        dN = 0.5*(dN + dN.transpose().eval());
        
        // Compute tangent-space basis from dN
        //   1 - pick the column with maximal norm as the first tangent vector,
        Index i0, i1, i2;
        Scalar sqNorm = dN.colwise().squaredNorm().maxCoeff(&i0);
        Mat32 B;
        B.col(0) = dN.col(i0) / sqrt(sqNorm);
        //   2 - orthogonalize the other column vectors, and pick the most reliable one
        i1 = (i0+1)%3;
        i2 = (i0+2)%3;
        VectorType v1 = dN.col(i1) - B.col(0).dot(dN.col(i1)) * B.col(0);
        VectorType v2 = dN.col(i2) - B.col(0).dot(dN.col(i2)) * B.col(0);
        Scalar v1norm2 = v1.squaredNorm();
        Scalar v2norm2 = v2.squaredNorm();
        if(v1norm2 > v2norm2) B.col(1) = v1 / sqrt(v1norm2);
        else                  B.col(1) = v2 / sqrt(v2norm2);
        
        // Compute the 2x2 matrix representing the shape operator by transforming dN to the basis B.
        // Recall that dN is a bilinear form, it thus transforms as follows:
        Mat22 S = B.transpose() * dN * B;
        
        Eigen::SelfAdjointEigenSolver<Mat22> eig2;
        eig2.computeDirect(S);
        
        if (eig2.info() != Eigen::Success){
          return UNDEFINED;
        }
        
        m_k1 = eig2.eigenvalues()(0);
        m_k2 = eig2.eigenvalues()(1);
        
        m_v1 = B * eig2.eigenvectors().col(0);
        m_v2 = B * eig2.eigenvectors().col(1);
        
        if(std::abs(m_k1)<std::abs(m_k2))
        {
          std::swap(m_k1, m_k2);
          std::swap(m_v1, m_v2);
        }
    }

    return bResult;
}
