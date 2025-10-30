/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _NFilter, typename T>
void
SphereFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();
    m_matA.setZero();
}

template < class DataPoint, class _NFilter, typename T>
bool
SphereFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                     const VectorType &localQ,
                                                     const DataPoint &attributes)
{
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        VectorA a;
#ifdef __CUDACC__
        a(0) = 1;
        a.template segment<DataPoint::Dim>(1) = localQ;
        a(DataPoint::Dim+1) = localQ.squaredNorm();
#else
        a << 1, localQ, localQ.squaredNorm();
#endif
        m_matA     += w * a * a.transpose();
        return true;
    }

    return false;
}


template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
SphereFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    // Compute status
    if(Base::finalize() != STABLE)
        return Base::m_eCurrentState;
    if(Base::getNumNeighbors() < DataPoint::Dim)
        return Base::m_eCurrentState = UNDEFINED;
    if (Base::algebraicSphere().isValid())
        Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
    else
        Base::m_eCurrentState = Base::getNumNeighbors() < 2*DataPoint::Dim ? UNSTABLE : STABLE;

    MatrixA matC;
    matC.setIdentity();
    matC.template topRightCorner<1,1>()    << -2;
    matC.template bottomLeftCorner<1,1>()  << -2;
    matC.template topLeftCorner<1,1>()     << 0;
    matC.template bottomRightCorner<1,1>() << 0;

    MatrixA invCpratt;
    invCpratt.setIdentity();
    invCpratt.template topRightCorner<1,1>()    << -0.5;
    invCpratt.template bottomLeftCorner<1,1>()  << -0.5;
    invCpratt.template topLeftCorner<1,1>()     << 0;
    invCpratt.template bottomRightCorner<1,1>() << 0;

    // Remarks:
    //   A and C are symmetric so all eigenvalues and eigenvectors are real
    //   we look for the minimal positive eigenvalue (eigenvalues may be negative)
    //   C^{-1}A is not symmetric
    //   calling Eigen::GeneralizedEigenSolver on (A,C) and Eigen::EigenSolver on C^{-1}A is equivalent
    //   C is not positive definite so Eigen::GeneralizedSelfAdjointEigenSolver cannot be used
#ifdef __CUDACC__
    m_solver.computeDirect(invCpratt * m_matA);
#else
    m_solver.compute(invCpratt * m_matA);
#endif
    VectorA eivals = m_solver.eigenvalues().real();
    int minId = -1;
    for(int i=0 ; i<DataPoint::Dim+2 ; ++i)
    {
    Scalar ev = eivals(i);
    if((ev>0) && (minId==-1 || ev<eivals(minId)))
    minId = i;
    }

    //mLambda = eivals(minId);
    VectorA vecU = m_solver.eigenvectors().col(minId).real();
    Base::m_uq = vecU[1+DataPoint::Dim];
    Base::m_ul = vecU.template segment<DataPoint::Dim>(1);
    Base::m_uc = vecU[0];

    Base::m_isNormalized = false;

    return Base::m_eCurrentState;
}
