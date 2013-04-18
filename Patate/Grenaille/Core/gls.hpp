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
