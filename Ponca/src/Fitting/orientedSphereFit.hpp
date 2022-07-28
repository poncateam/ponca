/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template < class DataPoint, class _WFunctor, typename T>
void
OrientedSphereFitImpl<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);

    // Setup fitting internal values
    m_sumDotPN = Scalar(0.0);
    m_sumDotPP = Scalar(0.0);
    m_nume     = Scalar(0.0);
    m_deno     = Scalar(0.0);
}

template<class DataPoint, class _WFunctor, typename T>
bool
OrientedSphereFitImpl<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                        const VectorType &localQ,
                                                        const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_sumDotPN += w * attributes.normal().dot(localQ);
        m_sumDotPP += w * localQ.squaredNorm();
        return true;
    }
    return false;
}


template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
OrientedSphereFitImpl<DataPoint, _WFunctor, T>::finalize ()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);
    PONCA_MULTIARCH_STD_MATH(max);
    PONCA_MULTIARCH_STD_MATH(abs);

    // Compute status
    if(Base::finalize() != STABLE  || Base::m_nbNeighbors < 3)
        return Base::m_eCurrentState;
    if (Base::algebraicSphere().isValid())
        Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
    else
        Base::m_eCurrentState = Base::m_nbNeighbors < 6 ? UNSTABLE : STABLE;

    // 1. finalize sphere fitting
    Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    Scalar invSumW = Scalar(1.)/Base::m_sumW;

    m_nume = (m_sumDotPN - invSumW * Base::m_sumP.dot(Base::m_sumN));
    Scalar den1 = invSumW * Base::m_sumP.dot(Base::m_sumP);
    m_deno = m_sumDotPP - den1;

    // Deal with degenerate cases
    if(abs(m_deno) < epsilon * max(m_sumDotPP, den1))
    {
        //plane
        Scalar s = Scalar(1.) / Base::m_ul.norm();
        Base::m_ul = s*Base::m_ul;
        Base::m_uc = s*Base::m_uc;
        Base::m_uq = Scalar(0.);
    }
    else
    {
        //Generic case
        Base::m_uq = Scalar(.5) * m_nume / m_deno;
        Base::m_ul = (Base::m_sumN - Base::m_sumP * (Scalar(2.) * Base::m_uq)) * invSumW;
        Base::m_uc = -invSumW * (Base::m_ul.dot(Base::m_sumP) + m_sumDotPP * Base::m_uq);
    }

    Base::m_isNormalized = false;

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
void
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);

    m_dSumN.setZero();

    m_dSumDotPN.setZero();
    m_dSumDotPP.setZero();
    m_dNume.setZero();
    m_dDeno.setZero();

    m_dUc.setZero();
    m_dUq.setZero();
    m_dUl.setZero();
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
bool
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::addLocalNeighbor(Scalar w,
                                                                           const VectorType &localQ,
                                                                           const DataPoint &attributes,
                                                                           ScalarArray &dw) {
    if( Base::addLocalNeighbor(w, localQ, attributes, dw) ) {

        m_dSumN     += attributes.normal() * dw;
        m_dSumDotPN += dw * attributes.normal().dot(localQ);
        m_dSumDotPP += dw * localQ.squaredNorm();

        return true;
    }
    return false;
}


template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::finalize()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    Base::finalize();
    // Test if base finalize end on a viable case (stable / unstable)
    if (this->isReady())
    {
        Scalar invSumW = Scalar(1.)/Base::m_sumW;

        Scalar nume  = Base::m_sumDotPN - invSumW*Base::m_sumP.dot(Base::m_sumN);
        Scalar deno  = Base::m_sumDotPP - invSumW*Base::m_sumP.dot(Base::m_sumP);

        // FIXME, the following product "Base::m_sumN.transpose() * m_dSumP" is prone to numerical cancellation
        // issues for spacial derivatives because, (d sum w_i P_i)/(d x) is supposed to be tangent to the surface whereas
        // "sum w_i N_i" is normal to the surface...
        m_dNume = m_dSumDotPN
            - invSumW*invSumW * ( Base::m_sumW * ( Base::m_sumN.transpose() * Base::m_dSumP + Base::m_sumP.transpose() * m_dSumN )
            - Base::m_dSumW*Base::m_sumP.dot(Base::m_sumN) );

        m_dDeno = m_dSumDotPP
            - invSumW*invSumW*(   Scalar(2.) * Base::m_sumW * Base::m_sumP.transpose() * Base::m_dSumP
            - Base::m_dSumW*Base::m_sumP.dot(Base::m_sumP) );

        m_dUq =  Scalar(.5) * (deno * m_dNume - m_dDeno * nume)/(deno*deno);

        // FIXME: this line is prone to numerical cancellation issues because dSumN and u_l*dSumW are close to each other.
        // If using two passes, one could directly compute sum( dw_i + (n_i - ul) ) to avoid this issue.
        m_dUl =  invSumW * ( m_dSumN - Base::m_ul*Base::m_dSumW - Scalar(2.)*(Base::m_dSumP * Base::m_uq + Base::m_sumP * m_dUq) );
        m_dUc = -invSumW*( Base::m_sumP.transpose() * m_dUl
            + Base::m_sumDotPP * m_dUq
            + Base::m_ul.transpose() * Base::m_dSumP
            + Base::m_uq * m_dSumDotPP
            + Base::m_dSumW * Base::m_uc);
    }

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename OrientedSphereDerImpl <DataPoint, _WFunctor, DiffType, T>::VectorArray
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::dNormal() const
{
  // Computes the derivatives of the normal of the sphere at the evaluation point.
  // Therefore, we must take into account the variation of the evaluation point when differentiating wrt space
  // i.e., normal(x) = grad / |grad|, with grad(x) = ul + 2 uq * x, and diff_x(grad) = dul + 2 uq I
  VectorArray dgrad = m_dUl;
  if(Base::isSpaceDer())
    dgrad.template rightCols<DataPoint::Dim>().diagonal().array() += Scalar(2)*Base::m_uq;
  Scalar norm  = Base::m_ul.norm();
  Scalar norm3 = norm*norm*norm;
  return dgrad / norm - Base::m_ul * (Base::m_ul.transpose() * dgrad) / norm3;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename OrientedSphereDerImpl <DataPoint, _WFunctor, DiffType, T>::ScalarArray
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::dPotential() const
{
  ScalarArray dfield = m_dUc;
  if(Base::isSpaceDer())
    dfield.template tail<DataPoint::Dim>() += Base::m_ul;
  return dfield;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
bool
OrientedSphereDerImpl<DataPoint, _WFunctor, DiffType, T>::applyPrattNorm()
{
    if(Base::isNormalized())
        return false; //need original parameters without Pratt Normalization

    PONCA_MULTIARCH_STD_MATH(sqrt);
    Scalar pn2    = Base::prattNorm2();
    Scalar pn     = sqrt(pn2);

    ScalarArray dpn2   = dprattNorm2();
    ScalarArray factor = Scalar(0.5) * dpn2 / pn;

    m_dUc = ( m_dUc * pn - Base::m_uc * factor ) / pn2;
    m_dUl = ( m_dUl * pn - Base::m_ul * factor ) / pn2;
    m_dUq = ( m_dUq * pn - Base::m_uq * factor ) / pn2;

    Base::applyPrattNorm();
    return true;
}
