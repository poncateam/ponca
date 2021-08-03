
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _WFunctor, typename T>
void
LeastSquareLine<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    // Setup primitive
    Base::resetPrimitive();
    Base::basisCenter() = _evalPos;

    // Setup fitting internal values
    m_sum         = Scalar(0.0);
    m_cog         = VectorType::Zero();
    m_cov         = MatrixType::Zero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
LeastSquareLine<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& _nei)
{
    VectorType q = _nei.pos() - Base::basisCenter();
    m_cog  +=  q;
    m_sum++;
    m_cov  +=  q * q.transpose();
    ++(Base::m_nbNeighbors);

    return true;  
}

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
LeastSquareLine<DataPoint, _WFunctor, T>::finalize ()
{
    /* handle specific configurations
     With less than 2 neighbors the fitting is undefined */
    if(m_sum == Scalar(0.) || Base::m_nbNeighbors < 2)
    {
      Base::resetPrimitive();
      Base::m_eCurrentState = UNDEFINED;
      return Base::m_eCurrentState;
    }

    // Finalize the centroid (still expressed in local basis)
    m_cog = m_cog/m_sum;

    // Center the covariance on the centroid
    m_cov = m_cov/m_sum - m_cog * m_cog.transpose();
  
    #ifdef __CUDACC__
        m_solver.computeDirect(m_cov);
    #else
        m_solver.compute(m_cov);
    #endif

    Base::m_eCurrentState = ( m_solver.info() == Eigen::Success ? STABLE : UNDEFINED );

    Base::setLine(m_cog, m_solver.eigenvectors().col(2).normalized());

    return Base::m_eCurrentState;
}
