/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint, class _NFilter, typename T>       // Outer template for the class
template <typename NF, std::enable_if_t<NF::isLocal, int>>   // Inner template of the method to enable project only if NF::isLocal
typename AlgebraicSphere<DataPoint, _NFilter, T>::VectorType // Return type
AlgebraicSphere<DataPoint, _NFilter, T>::project(const VectorType& _q) const
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    // turn to centered basis
    const VectorType lq = Base::m_nFilter.convertToLocalBasis(_q);

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
    return Base::m_nFilter.convertToGlobalBasis( lq + t * grad );
}

template < class DataPoint, class _NFilter, typename T>
typename AlgebraicSphere<DataPoint, _NFilter, T>::VectorType
AlgebraicSphere<DataPoint, _NFilter, T>::projectDescent( const VectorType& _q, int nbIter) const
{
    PONCA_MULTIARCH_STD_MATH(min)

    // turn to centered basis
    const VectorType lq = Base::m_nFilter.convertToLocalBasis(_q);

    VectorType grad;
    VectorType dir  = primitiveGradient(lq, false);
    Scalar ilg      = Scalar(1.)/dir.norm();
    dir             = dir*ilg;
    Scalar ad       = potential(lq, false);
    Scalar delta    = -ad*min(ilg,Scalar(1.));
    VectorType proj = lq + dir*delta;

    for (int i=0; i<nbIter; ++i)
    {
        grad  = primitiveGradient(proj, false);
        ilg   = Scalar(1.)/grad.norm();
        delta = -potential(proj, false)*min(ilg,Scalar(1.));
        proj += dir*delta;
    }
    return Base::m_nFilter.convertToGlobalBasis( proj );
}

template < class DataPoint, class _NFilter, typename T>
typename AlgebraicSphere<DataPoint, _NFilter, T>::Scalar
AlgebraicSphere<DataPoint, _NFilter, T>::potential( const VectorType &_q, const bool convertToLocalBasis ) const
{
    // turn to centered basis
    const VectorType lq = convertToLocalBasis? Base::m_nFilter.convertToLocalBasis(_q) : _q;
    return m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
}


template < class DataPoint, class _NFilter, typename T>
typename AlgebraicSphere<DataPoint, _NFilter, T>::VectorType
AlgebraicSphere<DataPoint, _NFilter, T>::primitiveGradient( const VectorType &_q, const bool convertToLocalBasis ) const
{
        // turn to centered basis
        const VectorType lq = convertToLocalBasis? Base::m_nFilter.convertToLocalBasis(_q) : _q;
        return (m_ul + Scalar(2.) * m_uq * lq);
}

