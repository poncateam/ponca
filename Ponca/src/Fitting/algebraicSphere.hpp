/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template <class DataPoint, class _NFilter, typename T>       // Outer template for the class
typename AlgebraicSphere<DataPoint, _NFilter, T>::VectorType // Return type
AlgebraicSphere<DataPoint, _NFilter, T>::project(const VectorType& _q) const
{
    PONCA_MULTIARCH_STD_MATH(sqrt);

    // turn to centered basis
    const VectorType lq = Base::getNeighborFilter().convertToLocalBasis(_q);



    if(isPlane())
    {
        Scalar sqnorm = m_ul.squaredNorm();
        return Base::getNeighborFilter().convertToGlobalBasis( lq - m_ul*(lq.dot(m_ul))/sqnorm );
    }
    else
    {
        Scalar potential = potentialLocal(lq);
        VectorType grad = primitiveGradientLocal(lq);
        Scalar norm = grad.norm();
        Scalar t = - (norm - sqrt(norm*norm - Scalar(4) * m_uq * potential)) / (Scalar(2) * m_uq * norm);
        return Base::getNeighborFilter().convertToGlobalBasis( lq + t * grad );
    }

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

