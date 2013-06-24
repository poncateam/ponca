/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/




template < class DataPoint, class _WFunctor, typename T>
void 
OrientedSphereFit<DataPoint, _WFunctor, T>::init(const VectorType& evalPos){
  
  // Setup primitive
  Base::resetPrimitive();
  Base::basisCenter() = evalPos;
    
  // Setup fitting internal values
  _sumP     = VectorType::Zero();
  _sumN     = VectorType::Zero();
  _sumDotPN = Scalar(0.0);
  _sumDotPP = Scalar(0.0);
  _sumW     = Scalar(0.0);
}

template < class DataPoint, class _WFunctor, typename T>
void 
OrientedSphereFit<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& nei){
    
  // centered basis
  VectorType q = nei.pos() - Base::basisCenter();
  
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
    Base::_uq = Scalar(0.);
  }else{
    invSumW = Scalar(1.)/_sumW;

    Base::_uq = Scalar(.5)* (_sumDotPN - invSumW*_sumP.dot(_sumN))
      / (_sumDotPP - invSumW*_sumP.dot(_sumP));
  }    

  Base::_ul = (_sumN-_sumP*(Scalar(2.)*Base::_uq))*invSumW;
  Base::_uc = -invSumW*(Base::_ul.dot(_sumP) + _sumDotPP*Base::_uq);
    
  // 2. Deal with plane case:
  if (fabs(Base::_uq)<Scalar(1e-9)){
    Scalar s = Scalar(1.) / Base::_ul.norm();
    Base::_ul = s*Base::_ul;
    Base::_uc = s*Base::_uc;
    Base::_uq = Scalar(0.);
  }

  Base::_isNormalized = false;
}




namespace internal{

  template < class DataPoint, class _WFunctor, typename T, int Type>
  void 
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::init(const VectorType& evalPos){
    Base::init(evalPos);

    _dSumN     = VectorArray::Zero();
    _dSumP     = VectorArray::Zero();
      
    _dSumDotPN = ScalarArray::Zero();
    _dSumDotPP = ScalarArray::Zero();
    _dSumW     = ScalarArray::Zero();
      
    _dUc       = ScalarArray::Zero();
    _dUq       = ScalarArray::Zero();
    _dUl       = VectorArray::Zero();
  }


  template < class DataPoint, class _WFunctor, typename T, int Type>
  void 
  OrientedSphereDer<DataPoint, _WFunctor, T, Type>::addNeighbor(const DataPoint  &nei){
    Base::addNeighbor(nei);

    int spaceId = (Type & FitScaleDer) ? 1 : 0;

    ScalarArray w;

    // centered basis
    VectorType q = nei.pos()-Base::basisCenter();

    // compute weight
    if (Type & FitScaleDer)
      w[0] = Base::_w.scaledw(q, nei);

    if (Type & FitSpaceDer){
      VectorType vw = Base::_w.spacedw(q, nei);
      for(unsigned int i = 0; i < DataPoint::Dim; i++)
	      w [spaceId+i] = vw[i];
    }

    // increment
    _dSumW     += w;
    _dSumP     += q * w;
    _dSumN     += nei.normal() * w;
    _dSumDotPN += w * nei.normal().dot(q);
    _dSumDotPP += w * q.squaredNorm();
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
      
      ScalarArray dNume = _dSumDotPN - invSumW*invSumW * ( Base::_sumW * (
			                                           Base::_sumN.transpose() * _dSumP +
			                                           Base::_sumP.transpose() * _dSumN ) 
			                                       - _dSumW*Base::_sumP.dot(Base::_sumN) );
      ScalarArray dDeno = _dSumDotPP - invSumW*invSumW*( Scalar(2.)*Base::_sumW * Base::_sumP.transpose()*_dSumP
			                  - _dSumW*Base::_sumP.dot(Base::_sumP) );

      _dUq =  Scalar(.5) * (deno * dNume - dDeno * nume)/(deno*deno);
      _dUl =  invSumW*((_dSumN - Scalar(2.)*(_dSumP*Base::_uq+Base::_sumP*_dUq)) - Base::_ul*_dSumW);
      _dUc = -invSumW*( Base::_sumP.transpose() * _dUl + 
                        Base::_sumDotPP * _dUq + Base::_ul.transpose() * _dSumP +
                        Base::_uq*_dSumDotPP + _dSumW*Base::_uc);
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
      
    ScalarArray dpn2   = dprattNorm2();
    ScalarArray factor = Scalar(0.5) * dpn2 / pn;	
    
    _dUc = ( _dUc * pn - Base::_uc * factor ) / pn2;
    _dUl = ( _dUl * pn - Base::_ul * factor ) / pn2;
    _dUq = ( _dUq * pn - Base::_uq * factor ) / pn2;
    
    Base::applyPrattNorm();
    return true;
  }
  
}// namespace internal
