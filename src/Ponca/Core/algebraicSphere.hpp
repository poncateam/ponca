/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
AlgebraicSphere<DataPoint, _WFunctor, T>::project(const VectorType& _q) const
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    // turn to centered basis
    const VectorType lq = _q-m_p;

    Scalar potential = m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
    VectorType grad = m_ul + Scalar(2) * m_uq * lq;
    Scalar norm = grad.norm();

    Scalar t;
    if(isPlane())
    {
        t = - potential / (norm*norm);
    }
    else
    {
        t = - (norm - sqrt(norm*norm - Scalar(4) * m_uq * potential)) / (Scalar(2) * m_uq * norm);
    }
    return m_p + lq + t * grad;
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
AlgebraicSphere<DataPoint, _WFunctor, T>::projectDescent( const VectorType& _q ) const
{
    PONCA_MULTIARCH_STD_MATH(min)

    // turn to centered basis
    const VectorType lq = _q-m_p;

    VectorType grad;
    VectorType dir  = m_ul+Scalar(2.)*m_uq*lq;
    Scalar ilg      = Scalar(1.)/dir.norm();
    dir             = dir*ilg;
    Scalar ad       = m_uc + m_ul.dot(lq) + m_uq * lq.squaredNorm();
    Scalar delta    = -ad*min(ilg,Scalar(1.));
    VectorType proj = lq + dir*delta;

    for (int i=0; i<16; ++i)
    {
        grad  = m_ul+Scalar(2.)*m_uq*proj;
        ilg   = Scalar(1.)/grad.norm();
        delta = -(m_uc + proj.dot(m_ul) + m_uq * proj.squaredNorm())*min(ilg,Scalar(1.));
        proj += dir*delta;
    }
    return proj + m_p;
}

template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::Scalar
AlgebraicSphere<DataPoint, _WFunctor, T>::potential( const VectorType &_q ) const
{
    // turn to centered basis
    const VectorType lq = _q-m_p;

    return m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
}


template < class DataPoint, class _WFunctor, typename T>
typename DataPoint::VectorType
AlgebraicSphere<DataPoint, _WFunctor, T>::primitiveGradient( const VectorType &_q ) const
{
        // turn to centered basis
        const VectorType lq = _q-m_p;
        return (m_ul + Scalar(2.f) * m_uq * lq);
}

