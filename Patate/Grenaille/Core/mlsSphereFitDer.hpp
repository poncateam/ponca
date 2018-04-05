

template < class DataPoint, class _WFunctor, typename T>
void
MlsSphereFitDer<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);

    m_d2Uc = Matrix::Zero(),
    m_d2Uq = Matrix::Zero();
    m_d2Ul = MatrixArray::Zero();

    m_d2SumDotPN = Matrix::Zero();
    m_d2SumDotPP = Matrix::Zero();
    m_d2SumW     = Matrix::Zero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
MlsSphereFitDer<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& _nei)
{
    bool bResult = Base::addNeighbor(_nei);

    if(bResult)
    {
        // centered basis
        VectorType q = _nei.pos() - Base::basisCenter();

        Matrix d2w;
        d2w;//TODO(thib) computes d2w wrt Scale/Space der type

        m_d2SumDotPN += d2w * _nei.normal().dot(q);
        m_d2SumDotPP += d2w * q.squaredNorm();
        m_d2SumW     += d2w;
    }

    return bResult;
}

//TODO(thib) delete this macro
template<typename T> struct ____The_Unkown_Type_Is____;
#define GLUE1(X,Y) X##Y
#define GLUE(X,Y) GLUE1(X,Y)
#define WhatIsTheTypeOf(expr) ____The_Unkown_Type_Is____<decltype(expr)> GLUE(_,__LINE__)

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
MlsSphereFitDer<DataPoint, _WFunctor, T>::finalize()
{
    Base::finalize();

    if (this->isReady())
    {
        //TODO(thib) use nume, deno, dNume, dDeno from base class to avoid calculating them twice

        Scalar invSumW = Scalar(1.)/Base::m_sumW;

        Scalar nume  = Base::m_sumDotPN - invSumW*Base::m_sumP.dot(Base::m_sumN);
        Scalar deno  = Base::m_sumDotPP - invSumW*Base::m_sumP.dot(Base::m_sumP);

        ScalarArray dNume = Base::m_dSumDotPN
            - invSumW*invSumW * ( Base::m_sumW * ( Base::m_sumN.transpose() * Base::m_dSumP + Base::m_sumP.transpose() * Base::m_dSumN )
            - Base::m_dSumW*Base::m_sumP.dot(Base::m_sumN) );

        ScalarArray dDeno = Base::m_dSumDotPP
            - invSumW*invSumW*(   Scalar(2.) * Base::m_sumW * Base::m_sumP.transpose() * Base::m_dSumP
            - Base::m_dSumW*Base::m_sumP.dot(Base::m_sumP) );

        Matrix d2Nume = m_d2SumDotPN
            - invSumW*invSumW*invSumW*invSumW*(
                    Base::m_sumW*Base::m_sumW*(/*Base::m_sumW*(0)*/   //TODO(thib) replace this
                                               /*+ */Base::m_dSumW.transpose()*(Base::m_sumN.transpose()*Base::m_dSumP + Base::m_sumP.transpose()*Base::m_dSumN)
                                               - (Base::m_sumP.transpose()*Base::m_sumN)*m_d2SumW.transpose()
                                               - (Base::m_dSumN.transpose()*Base::m_sumP + Base::m_dSumP.transpose()*Base::m_sumN)*Base::m_dSumW)
                    - Scalar(2.)*Base::m_sumW*Base::m_dSumW.transpose()*(Base::m_sumW*(Base::m_sumN.transpose()*Base::m_dSumP + Base::m_sumP.transpose()*Base::m_dSumN)
                                                                         - (Base::m_sumP.transpose()*Base::m_sumN)*Base::m_dSumW));

        Matrix d2Deno = m_d2SumDotPP
            - invSumW*invSumW*invSumW*invSumW*(
                Base::m_sumW*Base::m_sumW*(/*Scalar(2.)*Base::m_sumW*(0)*/ //TODO(thib) replace this
                                           /*+ */Scalar(2.)*Base::m_dSumW.transpose()*(Base::m_sumP.transpose()*Base::m_dSumP)
                                           - (Base::m_sumP.transpose()*Base::m_sumP)*m_d2SumW.transpose()
                                           - Scalar(2.)*(Base::m_dSumP.transpose()*Base::m_sumP)*Base::m_dSumW)
                - Scalar(2.)*Base::m_sumW*Base::m_dSumW.transpose()*(Scalar(2.)*Base::m_sumW*Base::m_sumP.transpose()*Base::m_dSumP
                                                                     - (Base::m_sumP.transpose()*Base::m_sumP)*Base::m_dSumW));

        Scalar deno2 = deno*deno;

        m_d2Uq = Scalar(.5)/(deno2*deno2)*(deno2*(dDeno.transpose()*dNume + deno*d2Nume
                                                  - dNume.transpose()*dDeno - nume*d2Deno)
                                           - Scalar(2.)*deno*dDeno.transpose()*(deno*dNume - nume*dDeno));

        m_d2Ul; //TODO(thib) compute this
        
        m_d2Uc = -invSumW*(/*0*/ // TODO(thib) replace this
            //+ 0              // TODO(thib) replace this
            //+ 0              // TODO(thib) replace this
            //+ 0              // TODO(thib) replace this
            /*+ */Base::m_dUq.transpose()*Base::m_dSumDotPP + Base::m_uq*m_d2SumDotPP
            + Base::m_dSumDotPP.transpose()*Base::m_dUq + m_d2Uq*Base::m_sumDotPP
            + Base::m_uc*m_d2SumW + Base::m_dUc.transpose()*Base::m_dSumW
            - Base::m_dSumW.transpose()*Base::m_dUc);
    }

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, typename T>
typename MlsSphereFitDer<DataPoint, _WFunctor, T>::ScalarArray
MlsSphereFitDer<DataPoint, _WFunctor, T>::dPotential() const
{
    //TODO(thib) handle spaceDer/scaleDer
    return Base::m_dUc + Base::m_ul;
}

template < class DataPoint, class _WFunctor, typename T>
typename MlsSphereFitDer<DataPoint, _WFunctor, T>::VectorArray
MlsSphereFitDer<DataPoint, _WFunctor, T>::dNormal() const
{
    //TODO(thib) handle spaceDer/scaleDer
    //TODO(thib) this is just the hessian for now, need to normalize by the norm and then differenciate
    return m_d2Uc + Base::m_dUl + Base::m_dUl.transpose() + Scalar(2.)*Base::m_uq*VectorArray::Ones();
}
