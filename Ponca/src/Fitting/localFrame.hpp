template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename LocalFrame<DataPoint, _WFunctor, T>::VectorType
LocalFrame<DataPoint, _WFunctor, T>::worldToLocalFrame (const VectorType& _q) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  if (ignoreTranslation)
    return B.transpose() * _q;
  else {
    // apply rotation and translation to get uv coordinates
    return B.transpose() * (Base::m_w.convertToLocalBasis(_q));
  }
}

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename LocalFrame<DataPoint, _WFunctor, T>::VectorType
LocalFrame<DataPoint, _WFunctor, T>::localFrameToWorld (const VectorType& _lq) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  if (ignoreTranslation)
    return B * _lq;
  else {
    return B * _lq + Base::m_w.basisCenter();
  }
}

template < class DataPoint, class _WFunctor, typename T>
void
LocalFrame<DataPoint, _WFunctor, T>::computeFrameFromNormalVector (const VectorType& _norm)
{
  // Creation of the vector 'a' non-collinear to the normal vector
  VectorType a;
  if (std::abs(_norm.x()) > std::abs(_norm.z())) {
      a = VectorType(-_norm.y(), _norm.x(), 0);
  } else {
      a = VectorType(0, -_norm.z(), _norm.y());
  }
  a.normalize();
  // Creation of the two vectors of the local frame (m_u and m_v) thanks to the cross product
  VectorType m_u = _norm.cross(a);
  VectorType m_v = _norm.cross(m_u);
  m_u.normalize();
  m_v.normalize();
  setFrameUV (m_u, m_v);
}
