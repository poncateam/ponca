
template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::w(const VectorType& q, 
					   const DataPoint&) const{
  Scalar d  = q.norm();  
  return (d <= _t) ? _wk.f(d/_t) : Scalar(0.);
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::VectorType
DistWeightFunc<DataPoint, WeightKernel>::spacedw(const VectorType& q, 
						 const DataPoint&) const{
  VectorType result = VectorType::Zero();
  if (q.norm() <= _t){
    for(unsigned int d = 0; d!= DataPoint::Dim; d++)
      result[d] = _wk.df(q[d]/_t) / _t;
  }

  return result;
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::scaledw(const VectorType& q, 
						 const DataPoint&) const{
  Scalar d  = q.norm();  
  return (d <= _t) ? ( - d*_wk.df(d/_t)/(_t*_t) ) : Scalar(0.);
}
