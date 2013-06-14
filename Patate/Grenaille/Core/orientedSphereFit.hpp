/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/




template < class DataPoint, class _WFunctor, typename T>
void 
OrientedSphereFit<DataPoint, _WFunctor, T>::init(const VectorType& evalPos){
  _p    = evalPos;

  _isNormalized = false;
    
  // initial values
  _sumP     = VectorType::Zero();
  _sumN     = VectorType::Zero();
  _sumDotPN = Scalar(0.0);
  _sumDotPP = Scalar(0.0);
  _sumW     = Scalar(0.0);
    
    
  _uc = Scalar(0.0);
  _ul = VectorType::Zero();
  _uq = Scalar(0.0);
}

template < class DataPoint, class _WFunctor, typename T>
void 
OrientedSphereFit<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& nei){
    
  // centered basis
  VectorType q = nei.pos()-_p;
  
  // compute weight
  Scalar w = _w.w(q, nei);  
  
  if (w > Scalar(0.)){
    
    // increment matrix
    _sumP     += q * w;
    _sumN     += nei.normal() * w;
    _sumDotPN += w * nei.normal().dot(q);
    _sumDotPP += w * q.squaredNorm();
    _sumW     += w;
  }
}


template < class DataPoint, class _WFunctor, typename T>
void
OrientedSphereFit<DataPoint, _WFunctor, T>::finalize (){
  MULTIARCH_STD_MATH(sqrt);

  // 1. finalize sphere fitting
  Scalar invSumW;
    
  if(_sumW == Scalar(0.)){ // handle planar configurations
    invSumW = 10000000;
    _uq = Scalar(0.);
  }else{
    invSumW = Scalar(1.)/_sumW;

    _uq = Scalar(.5)* (_sumDotPN - invSumW*_sumP.dot(_sumN))
      / (_sumDotPP - invSumW*_sumP.dot(_sumP));
  }
    

  _ul = (_sumN-_sumP*(Scalar(2.)*_uq))*invSumW;
  _uc = -invSumW*(_ul.dot(_sumP) + _sumDotPP*_uq);
    
  // 2. Deal with plane case:
  if (fabs(_uq)<Scalar(1e-9)){
    Scalar s = Scalar(1.) / _ul.norm();
    _ul = s*_ul;
    _uc = s*_uc;
    _uq = Scalar(0.);
  }

  _isNormalized = false;
}


/*!
*/
template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
OrientedSphereFit<DataPoint, _WFunctor, T>::project( VectorType q ) const{
  MULTIARCH_STD_MATH(min)

  // centered basis
  q = q-_p;

  //if(_isPlane)
  //{
  VectorType grad;
  VectorType dir  = _ul+Scalar(2.)*_uq*q;
  Scalar ilg      = Scalar(1.)/dir.norm();
  dir             = dir*ilg;
  Scalar ad       = _uc + _ul.dot(q) +
    _uq * q.squaredNorm();
  Scalar delta    = -ad*min(ilg,Scalar(1.));
  VectorType proj = q + dir*delta;

  for (int i=0 ; i<16 ; ++i)
    {
      grad = _ul+Scalar(2.)*_uq*proj;
      ilg = Scalar(1.)/grad.squaredNorm();
      delta = -(_uc + _ul.dot(proj) +
		_uq * proj.squaredNorm())*min(ilg,1.);
      proj += dir*delta;
    }
  return proj + _p;
  //}
  //return other - _ul * dot(other,_ul) + _uc;
  //return normalize(other-_center) * _r + _center;
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::Scalar
OrientedSphereFit<DataPoint, _WFunctor, T>::evaluate( VectorType q ) const{
  return _uc + q.dot(_ul) + _uq * q.dot(q);
}


namespace internal{

  template < class DataPoint, class _WFunctor, typename T, int Type>
  void 
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::init(const VectorType& evalPos){
    Base::init(evalPos);

    for (unsigned int d = 0; d < derDimension(); d++){

      _dSumN[d]     = VectorType::Zero();
      _dSumP[d]     = VectorType::Zero();
        
      _dSumDotPN[d] = Scalar(0.0);
      _dSumDotPP[d] = Scalar(0.0);
      _dSumW[d]     = Scalar(0.0);
        
      _dUc[d] = Scalar(0.0);
      _dUq[d] = Scalar(0.0);
      _dUl[d] = VectorType::Zero();
    }
  }


  template < class DataPoint, class _WFunctor, typename T, int Type>
  void 
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::addNeighbor(const DataPoint  &nei){
    Base::addNeighbor(nei);

    int spaceId = (Type & FitScaleDer) ? 1 : 0;

    ScalarArray w;

    // centered basis
    VectorType q = nei.pos()-Base::_p;

    // compute weight
    if (Type & FitScaleDer)
      w[0] = Base::_w.scaledw(q, nei);

    if (Type & FitSpaceDer){
      VectorType vw = Base::_w.spacedw(q, nei);
      for(unsigned int i = 0; i < DataPoint::Dim; i++)
	w [spaceId+i] = vw[i];
    }

    // increment
    for (unsigned int d = 0; d < derDimension(); d++){
      _dSumW[d]     += w[d];
      _dSumP[d]     += w[d] * q;
      _dSumN[d]     += w[d] * nei.normal();
      _dSumDotPN[d] += w[d] * nei.normal().dot(q);
      _dSumDotPP[d] += w[d] * q.squaredNorm();
    }
  }


  template < class DataPoint, class _WFunctor, typename T, int Type>
  void 
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::finalize(){
    MULTIARCH_STD_MATH(sqrt);

    Base::finalize();
    
    if (Base::_sumW != Scalar(0.)){

      Scalar invSumW = Scalar(1.)/Base::_sumW;

      Scalar nume  = Base::_sumDotPN - invSumW*Base::_sumP.dot(Base::_sumN);
      Scalar deno  = Base::_sumDotPP - invSumW*Base::_sumP.dot(Base::_sumP);

      // increment
      for (unsigned int d = 0; d < derDimension(); d++){
	Scalar dNume = _dSumDotPN[d] - invSumW*invSumW*(
							Base::_sumW*(_dSumP[d].dot(Base::_sumN)+Base::_sumP.dot(_dSumN[d])) - _dSumW[d]*Base::_sumP.dot(Base::_sumN));
	Scalar dDeno = _dSumDotPP[d] - invSumW*invSumW*( Scalar(2.)*Base::_sumW*_dSumP[d].dot(Base::_sumP)
							 - _dSumW[d]*Base::_sumP.dot(Base::_sumP));

	_dUq[d] = Scalar(.5) * (deno * dNume - dDeno * nume)/(deno*deno);
	_dUl[d] = invSumW*((_dSumN[d] - Scalar(2.)*(_dSumP[d]*Base::_uq+Base::_sumP*_dUq[d])) - _dSumW[d]*Base::_ul);
	_dUc[d] = -invSumW*( _dUl[d].dot(Base::_sumP) + _dUq[d]*Base::_sumDotPP + Base::_ul.dot(_dSumP[d])
			     + Base::_uq*_dSumDotPP[d] + _dSumW[d]*Base::_uc);
      }
    }

  }

  
  template < class DataPoint, class _WFunctor, typename T, int Type>
  bool
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::applyPrattNorm() {
    if(Base::isNormalized())
      return false; //need original parameters without Pratt Normalization

    MULTIARCH_STD_MATH(sqrt);
    Scalar pn2    = Base::prattNorm2();
    Scalar pn     = sqrt(pn2);
    
    for (unsigned int d = 0; d < derDimension(); d++){
      Scalar dpn2   = dprattNorm2(d);
      Scalar factor = Scalar(0.5) * dpn2 / pn;	
      
      _dUc[d] = ( _dUc[d] * pn - Base::_uc * factor ) / pn2;
      _dUl[d] = ( _dUl[d] * pn - Base::_ul * factor ) / pn2;
      _dUq[d] = ( _dUq[d] * pn - Base::_uq * factor ) / pn2;
    }
    
    Base::applyPrattNorm();
    return true;
  }
  
}// namespace internal
