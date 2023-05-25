

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename MeanPlane<DataPoint, _WFunctor, T>::VectorType
MeanPlane<DataPoint, _WFunctor, T>::worldToTangentPlane (const VectorType& _q) const
{
  if (ignoreTranslation)
    return B.transpose() * _q;
  else {
    // apply rotation and translation to get uv coordinates
    return B.transpose() * (Base::m_w.convertToLocalBasis(_q));
  }
}

template < class DataPoint, class _WFunctor, typename T>
template <bool ignoreTranslation>
typename MeanPlane<DataPoint, _WFunctor, T>::VectorType
MeanPlane<DataPoint, _WFunctor, T>::tangentPlaneToWorld (const VectorType& _lq) const
{
  if (ignoreTranslation)
    return B * _lq;
  else {
    return B * _lq + Base::m_w.basisCenter();
  }
}