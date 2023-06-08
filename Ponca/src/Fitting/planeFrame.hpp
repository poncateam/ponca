template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename PlaneFrame<DataPoint, _WFunctor, T>::VectorType
PlaneFrame<DataPoint, _WFunctor, T>::worldToLocalFrame (const VectorType& _q) const
{
  MatrixType B;
  B << Base::plane().normal(), m_u, m_v;
  if (ignoreTranslation)
    return B.transpose() * _q;
  else {
    // apply rotation and translation to get uv coordinates
    return B.transpose() * (Base::m_w.convertToLocalBasis(_q));
  }
}

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename PlaneFrame<DataPoint, _WFunctor, T>::VectorType
PlaneFrame<DataPoint, _WFunctor, T>::localFrameToWorld (const VectorType& _lq) const
{
  MatrixType B;
  B << Base::plane().normal(), m_u, m_v;
  if (ignoreTranslation)
    return B * _lq;
  else {
    return B * _lq + Base::m_w.basisCenter();
  }
}
