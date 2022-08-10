/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _WFunctor, typename T>
void
CovarianceFitBase<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);
    m_cov.setZero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
CovarianceFitBase<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                              const VectorType &localQ,
                                                              const DataPoint &attributes)
{
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_cov  += localQ * localQ.transpose();
        return true;
    }
    return false;
}


template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
CovarianceFitBase<DataPoint, _WFunctor, T>::finalize ()
{
    // handle specific configurations
    // With less than 3 neighbors the fitting is undefined
    if(Base::finalize() != STABLE || Base::m_nbNeighbors < 3)
    {
        return Base::m_eCurrentState = UNDEFINED;
    }

    // Center the covariance on the centroid
    auto centroid = Base::barycenter();
    m_cov = m_cov/Base::m_sumW - centroid * centroid.transpose();

#ifdef __CUDACC__
    m_solver.computeDirect(m_cov);
#else
    m_solver.compute(m_cov);
#endif
    Base::m_eCurrentState = ( m_solver.info() == Eigen::Success ? STABLE : UNDEFINED );

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, typename T>
typename CovarianceFitBase<DataPoint, _WFunctor, T>::Scalar
CovarianceFitBase<DataPoint, _WFunctor, T>::surfaceVariation () const
{
    return m_solver.eigenvalues()(0) / m_solver.eigenvalues().mean();
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
void
CovarianceFitDer<DataPoint, _WFunctor, DiffType, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);

    for(int k=0; k<Base::NbDerivatives; ++k)
        m_dCov[k].setZero();
}



template < class DataPoint, class _WFunctor, int DiffType, typename T>
bool
CovarianceFitDer<DataPoint, _WFunctor, DiffType, T>::addLocalNeighbor(Scalar w,
                                                                      const VectorType &localQ,
                                                                      const DataPoint &attributes,
                                                                      ScalarArray &dw)
{
    if( Base::addLocalNeighbor(w, localQ, attributes, dw) ) {
        for(int k=0; k<Base::NbDerivatives; ++k)
            m_dCov[k]  += dw[k] * localQ * localQ.transpose(); /// \fixme better use eigen here

        return true;
    }

    return false;
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
CovarianceFitDer<DataPoint, _WFunctor, DiffType, T>::finalize()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    Base::finalize();
    // Test if base finalize end on a viable case (stable / unstable)
    if (this->isReady())
    {
        VectorType cog = Base::barycenter();
        MatrixType cogSq = cog * cog.transpose();
        // \fixme Better use eigen here
        for(int k=0; k<Base::NbDerivatives; ++k)
        {
            // Finalize the computation of dCov.
            m_dCov[k] = m_dCov[k]
                        - cog * Base::m_dSumP.col(k).transpose()
                        - Base::m_dSumP.col(k) * cog.transpose()
                        + Base::m_dSumW[k] * cogSq;
        }
    }

    return Base::m_eCurrentState;
}

