template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau() const{
  MULTIARCH_STD_MATH(sqrt);
  
  ScalarArray result;
  Scalar prattNorm2 = Base::prattNorm2();
  Scalar prattNorm  = sqrt(prattNorm2);
  Scalar cfactor    = Scalar(.5) / prattNorm;
  
  for (unsigned int d = 0; d < Base::derDimension(); d++)
    result[d] = (Base::_dUc[d] * prattNorm - Base::_uc * cfactor * Base::dprattNorm2(d)) / prattNorm2;
  return result;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta() const{
  MULTIARCH_STD_MATH(sqrt);
  
  VectorArray result;
  Scalar prattNorm2 = Base::prattNorm2();
  Scalar prattNorm  = sqrt(prattNorm2);
  Scalar cfactor    = Scalar(.5) / prattNorm;
  Scalar deno       = Scalar(1.) / prattNorm2;
  
  for (unsigned int d = 0; d < Base::derDimension(); d++)
    result[d] = (Base::_dUl[d] * prattNorm - Base::_ul * cfactor * Base::dprattNorm2(d) ) * deno;
  return result;

}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa() const{
  MULTIARCH_STD_MATH(sqrt);
  
  ScalarArray result;
  Scalar prattNorm2 = Base::prattNorm2();
  Scalar prattNorm  = sqrt(prattNorm2);
  Scalar cfactor    = Scalar(.5) / prattNorm;
  
  for (unsigned int d = 0; d < Base::derDimension(); d++)
    result[d] = Scalar(2.) * (Base::_dUq[d] * prattNorm - Base::_uq * cfactor * Base::dprattNorm2(d)) / prattNorm2;
  return result;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau_normalized() const{
  return dtau();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta_normalized() const{
  Scalar ulN  = Base::_ul.norm();
  Scalar ulN3 = ulN*ulN*ulN;
  VectorArray detaStatic = deta();
  VectorArray result;
  for (unsigned int d = 0; d < Base::derDimension(); d++)
    result[d] = Base::_t * detaStatic[d]*(Scalar(1.)/ulN) - Base::_ul * detaStatic[d].dot(Base::_ul)*(Scalar(1.)/ulN3);
  return result;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa_normalized() const{
  return dkappa()*Base::_t*Base::_t;
}




template < class DataPoint, class _WFunctor, typename T>
typename GLSGeomVar <DataPoint, _WFunctor, T>::Scalar
GLSGeomVar <DataPoint, _WFunctor, T>::geomVar(Scalar wtau, 
							     Scalar weta,
							     Scalar wkappa) const{
  Scalar dtau   = Base::dtau_normalized()[0];
  Scalar deta   = Base::deta_normalized()[0].norm();
  Scalar dkappa = Base::dkappa_normalized()[0];
  return wtau*dtau*dtau + weta*deta*deta + wkappa*dkappa*dkappa;
  
}
