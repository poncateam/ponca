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
