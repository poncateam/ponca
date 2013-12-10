/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/



#if neverdefined
    MatrixBB cov = MatrixBB::Zero();
    VectorDd sumP = VectorDd::Zero();
    double sumDotPP = 0.;
    double sumOfWeights = 0.;
    
    // the normalization matrix
    MatrixBB Q = MatrixBB::Zero();
    
    for(uint i=0 ; i<nofSamples ; ++i)
    {
      VectorDd p = pNeighborhood->getNeighbor(i).position().cast<double>();
      VectorDd n = pNeighborhood->getNeighbor(i).normal().cast<double>();
      double w = pNeighborhood->getNeighborWeight(i);

      VectorB basis;
      basis << n, n.dot(p);
      
      cov += w * basis * basis.transpose();
      sumOfWeights += w;
      
      MatrixBB q;
      q <<  MatrixDDd::Identity(), p,
            p.transpose(), p.squaredNorm();
      
      Q += w * q;
      
      sumP      += w * p;
      sumDotPP  += w * p.squaredNorm();
    }
    cov /= sumOfWeights;
    Q   /= sumOfWeights;            
    
    MatrixBB M = Q.inverse() * cov;
    Eigen::EigenSolver<MatrixBB> eig(M);
    VectorB eivals = eig.eigenvalues().real();
    int maxId = 0;
    double l = eivals.maxCoeff(&maxId);
    VectorB eivec = eig.eigenvectors().col(maxId).real();
    
    // integrate
    uLinear()   = eivec.start<Dim>().cast<Real>();
    uQuad()     = 0.5*eivec(Dim);
    uConstant() = -(1./sumOfWeights)*(eivec.start<Dim>().dot(sumP) + 0.5*eivec(Dim) * sumDotPP);
#endif

template < class DataPoint, class _WFunctor, typename T>
void 
UnorientedSphereFit<DataPoint, _WFunctor, T>::init(const VectorType& evalPos){
  
  // Setup primitive
  Base::resetPrimitive();
  Base::basisCenter() = evalPos;
    
  // Setup fitting internal values
  _matA.setZero();
//   _matQ.setZero();
  _sumP.setZero();
  _sumDotPP = Scalar(0.0);
  _sumW     = Scalar(0.0);
}

template < class DataPoint, class _WFunctor, typename T>
void 
UnorientedSphereFit<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& nei){
    
  // centered basis
  VectorType q = nei.pos() - Base::basisCenter();
  
  // compute weight
  Scalar w = _w.w(q, nei);  
  
  if (w > Scalar(0.)){
    
    VectorB basis;
    basis << nei.normal(), nei.normal().dot(q);
    
    _matA     += w * basis * basis.transpose();
    _sumP     += w * q;
    _sumDotPP += w * q.squaredNorm();
    _sumW     += w;
  }
}


template < class DataPoint, class _WFunctor, typename T>
void
UnorientedSphereFit<DataPoint, _WFunctor, T>::finalize (){
  MULTIARCH_STD_MATH(sqrt);

  // 1. finalize sphere fitting
  Scalar invSumW;
  Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();
    
  if(_sumW == Scalar(0.)){ // handle empty configurations
    Base::_ul.setZero();
    Base::_uc = 0;
    Base::_uq = 0;
    Base::_isNormalized = false;
	Base::_isReady = false;
    return;
  }else{
    invSumW = Scalar(1.)/_sumW;
  }
  
  MatrixBB Q;
  Q.template topLeftCorner<Dim,Dim>().setIdentity();
  Q.col(Dim).template head<Dim>() = _sumP*invSumW;
  Q.row(Dim).template head<Dim>() = _sumP*invSumW;
  Q(Dim,Dim) = _sumDotPP*invSumW;
  _matA *= invSumW;
  
  MatrixBB M = Q.inverse() * _matA;
  Eigen::EigenSolver<MatrixBB> eig(M);
  VectorB eivals = eig.eigenvalues().real();
  int maxId = 0;
  Scalar l = eivals.maxCoeff(&maxId);
  VectorB eivec = eig.eigenvectors().col(maxId).real();
  
  // integrate
  Base::_ul = eivec.template head<Dim>();
  Base::_uq = 0.5*eivec(Dim);
  Base::_uc = -invSumW*(Base::_ul.dot(_sumP) + _sumDotPP*Base::_uq);
    
  Base::_isNormalized = false;
  Base::_isReady      = true;
}

#ifdef TOBEIMPLEMENTED

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

#endif
