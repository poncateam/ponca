/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>


 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template < class DataPoint, class _WFunctor, typename T>
void
LeastSquareSurfaceFit<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    // Setup primitive
    Base::resetPrimitive();
    Base::basisCenter() = _evalPos;

    // Setup fitting internal values
    m_sumW        = Scalar(0.0);
    m_cog         = VectorType::Zero();
    cofficient    = Eigen::Matrix<Scalar, 9, 1> ::Zero();
    m_right       = Eigen::Matrix<Scalar, 9, 1> ::Zero();
    m_cov         = Eigen::Matrix<Scalar, 9, 9> ::Zero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
LeastSquareSurfaceFit<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& _nei)
{
    VectorType q = _nei.pos() - Base::basisCenter();

    // compute weight
    Scalar w = m_w.w(q, _nei);

    if (w > Scalar(0.))
    {
      
      m_sumW  += w;
      m_cog  +=  q * w;

      // temp = [x^2 y^2 z^2 xy yz zx x y z]
      Eigen::Matrix<Scalar, 9, 1> temp;
      temp << q[0] * q[0], q[1] * q[1] ,q[2] * q[2] , 2* q[0] * q[1] , 2*q[1] * q[2] , 2 * q[0] * q[2] , 2*q; 
     
      m_right -= temp;
      m_cov  +=  temp * temp.transpose(); 

      ++(Base::m_nbNeighbors);
      
      return true;
    }

    return false;
}

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
LeastSquareSurfaceFit<DataPoint, _WFunctor, T>::finalize ()
{

    // handle specific configurations
    // With less than 2 neighbors the fitting is undefined
    if(m_sumW == Scalar(0.) || Base::m_nbNeighbors < 2)
    {
      Base::resetPrimitive();
      Base::m_eCurrentState = UNDEFINED;
      return Base::m_eCurrentState;
    }

    m_cog = m_cog/m_sumW; 
   
    #ifdef __CUDACC__
      m_solver.computeDirect(m_cov);
    #else
      cofficient = m_cov.colPivHouseholderQr().solve(m_right) ;
    #endif

    Base::m_eCurrentState = ( m_right.isApprox(m_cov * cofficient)  == true ? STABLE : UNDEFINED );
    Base::setSurface(cofficient , m_cog);  // RETURN(CALC_QUADRIC_E"QUATION(LAMBDA,EIGEN_VECTORS))

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::Vectortype
LeastSquareSurfaceFit<DataPoint, _WFunctor, T>::primitiveGradient( const VectorType &_q ) const
{
      MatrixType A;
      A << cofficient[0], cofficient[5], cofficient[4],
           cofficient[5], cofficient[1], cofficient[3],
           cofficient[4], cofficient[3], cofficient[2];
      
      VectorType B;
      B << cofficient[6], cofficient[7], cofficient[8];
       VectorType result = A * _q + B;
      return result;

}