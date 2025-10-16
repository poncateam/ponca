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
    const VectorType lq = Base::getNeighborFilter().convertToLocalBasis(_q);

    Scalar potential = potentialLocal(lq);
    VectorType grad = primitiveGradientLocal(lq);
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
    return Base::getNeighborFilter().convertToGlobalBasis( lq + t * grad );
}

template < class DataPoint, class _NFilter, typename T>
typename AlgebraicSphere<DataPoint, _NFilter, T>::Scalar
AlgebraicSphere<DataPoint, _NFilter, T>::potentialLocal( const VectorType &_lq ) const
{
    return m_uc + _lq.dot(m_ul) + m_uq * _lq.squaredNorm();
}

template < class DataPoint, class _NFilter, typename T>
typename AlgebraicSphere<DataPoint, _NFilter, T>::VectorType
AlgebraicSphere<DataPoint, _NFilter, T>::primitiveGradientLocal( const VectorType &_lq ) const
{
        return (m_ul + Scalar(2.) * m_uq * _lq);
}

