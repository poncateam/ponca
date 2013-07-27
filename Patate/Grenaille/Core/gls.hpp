/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau() const{
  MULTIARCH_STD_MATH(sqrt);
  
  Scalar prattNorm2 = Base::prattNorm2();
  Scalar prattNorm  = sqrt(prattNorm2);
  Scalar cfactor    = Scalar(.5) / prattNorm;
  
  return (Base::_dUc * prattNorm - Base::_uc * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta() const{
  Scalar ulN  = Base::_ul.norm();
  Scalar ulN3 = ulN*ulN*ulN;
  return Base::_dUl*(Scalar(1.)/ulN) - Base::_ul * Base::_ul.transpose()*Base::_dUl*(Scalar(1.)/ulN3);
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa() const{
  MULTIARCH_STD_MATH(sqrt);
  
  Scalar prattNorm2 = Base::prattNorm2();
  Scalar prattNorm  = sqrt(prattNorm2);
  Scalar cfactor    = Scalar(.5) / prattNorm;
  
  return Scalar(2.) * (Base::_dUq * prattNorm - Base::_uq * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau_normalized() const{
  return dtau();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta_normalized() const{
  return Base::_t * deta();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa_normalized() const{
  return dkappa()*Base::_t*Base::_t;
}




template < class DataPoint, class _WFunctor, typename T>
typename GLSGeomVar <DataPoint, _WFunctor, T>::Scalar
GLSGeomVar <DataPoint, _WFunctor, T>::geomVar( Scalar wtau, 
							                                 Scalar weta,
							                                 Scalar wkappa ) const{
  Scalar dtau   = Base::dtau_normalized().col(0)(0);
  Scalar deta   = Base::deta_normalized().col(0).norm();
  Scalar dkappa = Base::dkappa_normalized().col(0)(0);
  
  return wtau*dtau*dtau + weta*deta*deta + wkappa*dkappa*dkappa;  
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSSpatialVariation<DataPoint, _WFunctor, T>::GLSSpatialEigenDecomposition
GLSSpatialVariation <DataPoint, _WFunctor, T>::projectedVariationDecomposition( 
                                                   Scalar wtau, 
		                                               Scalar weta,
                      			                       Scalar wkappa ) const{
  MULTIARCH_STD_MATH(acos);                    			        
  
  // rotation matrix to express vectors in tangent plane
  MatrixType localBasis;  
  VectorType crossVector = VectorType::Ones();
  crossVector.normalize();
  
  localBasis.col(2) = Base::eta();
  localBasis.col(1) = localBasis.col(2).cross(crossVector).normalized();
  localBasis.col(0) = localBasis.col(2).cross(localBasis.col(1) ).normalized();  
  
  // Init result structure, jacobian matrix and solver
  GLSSpatialEigenDecomposition result;    
  Eigen::Matrix<Scalar, int(DataPoint::Dim)+2, int(DataPoint::Dim)> jacobian;  
  Eigen::Matrix<Scalar, int(DataPoint::Dim)+2, int(DataPoint::Dim)-1> projJacobian;
  
  Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Scalar, 
                                               int(DataPoint::Dim)-1, 
                                               int(DataPoint::Dim)-1> > eigensolver;
  
  const unsigned int firstId = Base::isScaleDer() ? 1 : 0;
  
  // Set jacobian matrix with GLS derivatives
  jacobian.template block< 1, int(DataPoint::Dim) >(0,0) = wtau *
  Base::dtau_normalized().template block< 1, int(DataPoint::Dim) >(firstId,0);
  
  jacobian.template block< int(DataPoint::Dim), int(DataPoint::Dim) >(1,0) = weta *
  Base::deta_normalized().template block< int(DataPoint::Dim), int(DataPoint::Dim) >(firstId,0);
  
  jacobian.template block< 1, int(DataPoint::Dim) >(int(DataPoint::Dim)+1,0) = wkappa *
  Base::dkappa_normalized().template block< 1, int(DataPoint::Dim) >(firstId,0);
  
  // project on tangeant plane
  projJacobian = (jacobian * localBasis).template block<int(DataPoint::Dim)+2,int(DataPoint::Dim)-1>(0,0);

  // compute tensor and extract eigen values/vectors
  eigensolver.compute(projJacobian.transpose() * projJacobian);
  
  if (eigensolver.info() == Eigen::Success){
    result.first  = eigensolver.eigenvalues();
    
    // get eigenvectors back and apply reverse rotation (from local to global)
    result.second = Eigen::Matrix<Scalar, DataPoint::Dim, DataPoint::Dim-1>::Zero();
    result.second.template block<int(DataPoint::Dim)-1,
                                 int(DataPoint::Dim)-1>(0,0) = eigensolver.eigenvectors();                  
    result.second = (result.second.transpose() * localBasis.transpose()).transpose();
  }
  
  return result;                      			                       
}
