template < class DataPoint, class _NFilter, typename T>
typename LocalFrame<DataPoint, _NFilter, T>::VectorType
LocalFrame<DataPoint, _NFilter, T>::worldToFrameLocal (const VectorType& _q) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  return B.transpose() * _q;
}

template < class DataPoint, class _NFilter, typename T>
typename LocalFrame<DataPoint, _NFilter, T>::VectorType
LocalFrame<DataPoint, _NFilter, T>::frameToWorldLocal (const VectorType& _lq) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  return B * _lq;
}
