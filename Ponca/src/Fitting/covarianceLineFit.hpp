
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _WFunctor, typename T>
void
CovarianceLineFit<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);
    m_cog         = VectorType::Zero();
    m_cov         = MatrixType::Zero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
CovarianceLineFit<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                        const VectorType &localQ,
                                                        const DataPoint &attributes)
{
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_cog  += w * localQ; /// \fixme Replace by MeanPosition
        m_cov  +=  localQ * localQ.transpose();
        return true;
    }
    return false;
}

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
CovarianceLineFit<DataPoint, _WFunctor, T>::finalize ()
{
    /* handle specific configurations
     With less than 2 neighbors the fitting is undefined */
    if(Base::finalize() != STABLE || Base::m_nbNeighbors < 2)
    {
      return Base::m_eCurrentState = UNDEFINED;
    }

    // Finalize the centroid (still expressed in local basis)
    m_cog = m_cog/Base::m_sumW;

    // Center the covariance on the centroid
    m_cov = m_cov/Base::m_sumW - m_cog * m_cog.transpose();
  
    #ifdef __CUDACC__
        m_solver.computeDirect(m_cov);
    #else
        m_solver.compute(m_cov);
    #endif

    Base::m_eCurrentState = ( m_solver.info() == Eigen::Success ? STABLE : UNDEFINED );

    Base::setLine(m_cog, m_solver.eigenvectors().col(2).normalized());

    return Base::m_eCurrentState;
}
