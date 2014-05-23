/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


template < class DataPoint, class _WFunctor, typename T>
void 
CovariancePlaneFit<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    // Setup primitive
    Base::resetPrimitive();
    // Base::basisCenter() = _evalPos;

    // Setup fitting internal values
    m_sumW        = Scalar(0.0);
    m_evalPos     = _evalPos;
    m_gc          = VectorType::Zero();
    m_cov         = MatrixType::Zero();
    m_isFirstPass = true;
}

template < class DataPoint, class _WFunctor, typename T>
bool 
CovariancePlaneFit<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& _nei)
{
    // compute weight
    Scalar w = m_w.w(_nei.pos() - m_evalPos, _nei);

    if (w > Scalar(0.))
    {
        if (m_isFirstPass) // compute gravity center
        { 
            m_gc   += w * _nei.pos();
            m_sumW += w;
        }
        else // increment covariance matrix
        {
            VectorType q = _nei.pos() - m_gc;
            m_cov +=  w * q * q.transpose();
        }
        
        ++(Base::m_nbNeighbors);
        return true;
    }
    return false;
}


template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
CovariancePlaneFit<DataPoint, _WFunctor, T>::finalize ()
{
    Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    // handle specific configurations
    // With less than 3 neighbors the fitting is undefined
    if(m_sumW == Scalar(0.) || Base::m_nbNeighbors < 3)
    {
        Base::resetPrimitive();
        Base::m_eCurrentState = UNDEFINED;
        return Base::m_eCurrentState;
    }

    if( m_isFirstPass )
    {
        m_gc  /= m_sumW;
        m_cov  = MatrixType::Zero();
        m_sumW = Scalar(0.);
        m_isFirstPass = false;
        
        Base::m_eCurrentState = NEED_OTHER_PASS;
        Base::m_nbNeighbors   = 0;
    }
    else // second pass
    {
        m_cov /= m_sumW;
        
#ifdef __CUDACC__
        m_solver.computeDirect(m_cov.transpose()*m_cov);
#else
        m_solver.compute(m_cov.transpose()*m_cov);
#endif
        
        Base::setPlane(m_solver.eigenvectors().col(0), m_gc);

        // \todo Use the output of the solver to check stability
        Base::m_eCurrentState = STABLE;
    }

    return Base::m_eCurrentState;
}


template < class DataPoint, class _WFunctor, typename T>
typename CovariancePlaneFit<DataPoint, _WFunctor, T>::Scalar
CovariancePlaneFit<DataPoint, _WFunctor, T>::surfaceVariation () const
{
    if( Base::m_eCurrentState = UNDEFINED )
      return 0;

    return m_solver.eigenvalues()(0) / m_solver.eigenvalues().norm();
}
