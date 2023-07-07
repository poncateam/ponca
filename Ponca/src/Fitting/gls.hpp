/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/

template < class DataPoint, class _WFunctor, typename T>
inline FIT_RESULT
GLSParam<DataPoint, _WFunctor, T>::finalize()
{
    FIT_RESULT bResult = Base::finalize();

    if(bResult != UNDEFINED)
    {
        m_fitness = Scalar(1.) - Base::prattNorm2();
    }

    return bResult;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::dtau() const
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;
    ScalarArray dfield = Base::m_dUc;
    // Recall that tau is the field function at the evaluation point, we thus must take care about
    // its variation when differentiating in space:
    if(this->isSpaceDer())
      dfield.template tail<DataPoint::Dim>() += Base::m_ul;

    return (dfield * prattNorm - Base::m_uc * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::VectorArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::deta() const
{
  return Base::dNormal();
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::dkappa() const
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;

    return Scalar(2.) * (Base::m_dUq * prattNorm - Base::m_uq * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::dtau_normalized() const
{
    return dtau();
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::VectorArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::deta_normalized() const
{
    return Base::m_w.evalScale() * deta();
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, DiffType, T>::dkappa_normalized() const
{
    return dkappa() * Base::m_w.evalScale() * Base::m_w.evalScale();
}




template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename GLSDer <DataPoint, _WFunctor, DiffType, T>::Scalar
GLSDer <DataPoint, _WFunctor, DiffType, T>::geomVar(  Scalar wtau,
                                                      Scalar weta,
                                                      Scalar wkappa ) const
{
    static_assert( bool(DiffType & FitScaleDer), "Scale derivatives are required to compute Geometric Variation" );
    Scalar dtau   = dtau_normalized().col(0)(0);
    Scalar deta   = deta_normalized().col(0).norm();
    Scalar dkappa = dkappa_normalized().col(0)(0);

    return wtau*dtau*dtau + weta*deta*deta + wkappa*dkappa*dkappa;
}
