/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _WFunctor, typename T>
void
SphereFitImpl<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);
    m_matA.setZero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
SphereFitImpl<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
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


template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
SphereFitImpl<DataPoint, _WFunctor, T>::finalize ()
{
    // Compute status
    if(Base::finalize() != STABLE || Base::m_nbNeighbors < 3)
        return Base::m_eCurrentState = UNDEFINED;
    if (Base::algebraicSphere().isValid())
        Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
    else
        Base::m_eCurrentState = Base::m_nbNeighbors < 6 ? UNSTABLE : STABLE;

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

    MatrixA M = invCpratt * m_matA;
    // go to positive semi-definite matrix to be compatible with
    // SelfAdjointEigenSolver requirements
    // Note: This does not affect the eigen vectors order
    Eigen::SelfAdjointEigenSolver<MatrixA> solver;
#ifdef __CUDACC__
    solver.computeDirect(M.transpose() * M);
#else
    solver.compute(M.transpose() * M);
#endif
    VectorA eivals = solver.eigenvalues().real();
    int minId = -1;
    for(int i=0 ; i<DataPoint::Dim+2 ; ++i)
    {
    Scalar ev = eivals(i);
    if((ev>0) && (minId==-1 || ev<eivals(minId)))
    minId = i;
    }

    //mLambda = eivals(minId);
    VectorA vecU = solver.eigenvectors().col(minId).real();
    Base::m_uq = vecU[1+DataPoint::Dim];
    Base::m_ul = vecU.template segment<DataPoint::Dim>(1);
    Base::m_uc = vecU[0];

    Base::m_isNormalized = false;

    return Base::m_eCurrentState;
}
