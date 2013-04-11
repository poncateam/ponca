
template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::w(const VectorType& relativeQuery, 
					   const DataPoint&){
  Scalar d  = relativeQuery.norm();  
  return (d <= _t) ? _wk.f(d/_t) : Scalar(0.);
}


template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::spacedw(const VectorType& relativeQuery, 
					   const DataPoint&){
  return Scalar(1.);
}


template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::scaledw(const VectorType& relativeQuery, 
					   const DataPoint&){
  return Scalar(1.);

}
