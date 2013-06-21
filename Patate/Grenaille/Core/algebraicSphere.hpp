
/*!
  Use gradient descent
*/
template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
AlgebraicSphere<DataPoint, _WFunctor, T>::project( const VectorType& q ) const{
  MULTIARCH_STD_MATH(min)

  // turn to centered basis
  const VectorType lq = q-_p;

  //if(_isPlane)
  //{
  VectorType grad;
  VectorType dir  = _ul+Scalar(2.)*_uq*lq;
  Scalar ilg      = Scalar(1.)/dir.norm();
  dir             = dir*ilg;
  Scalar ad       = _uc + _ul.dot(lq) + _uq * lq.squaredNorm();
  Scalar delta    = -ad*min(ilg,Scalar(1.));
  VectorType proj = lq + dir*delta;

  for (int i=0 ; i<16 ; ++i)
    {
      grad  = _ul+Scalar(2.)*_uq*proj;
      ilg   = Scalar(1.)/grad.norm();
      delta = -(_uc + proj.dot(_ul) + _uq * proj.squaredNorm())*min(ilg,Scalar(1.));
      proj += dir*delta;
    }
  return proj + _p;
  //}
  //return other - _ul * dot(other,_ul) + _uc;
  //return normalize(other-_center) * _r + _center;
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::Scalar
AlgebraicSphere<DataPoint, _WFunctor, T>::potential( const VectorType &q ) const{  
  // turn to centered basis
  const VectorType lq = q-_p;
  
  return _uc + lq.dot(_ul) + _uq * lq.squaredNorm();
}


template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
AlgebraicSphere<DataPoint, _WFunctor, T>::primitiveGradient( const VectorType &q ) const{
  // turn to centered basis
  const VectorType lq = q-_p;  
  return (_ul + Scalar(2.f) * _uq * lq);
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::MatrixType
AlgebraicSphere<DataPoint, _WFunctor, T>::primitiveHessian( const VectorType &q ) const{
	// this is very approximate !!
	return Scalar(2.) * _uq * MatrixType::Identity();	
}


