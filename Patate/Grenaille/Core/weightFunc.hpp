/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/



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
  Scalar d = q.norm();
  if (d <= _t) result = (q / (d * _t)) * _wk.df(d/_t);
  return result;
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::scaledw(const VectorType& q, 
						 const DataPoint&) const{
  Scalar d  = q.norm();  
  return (d <= _t) ? ( - d*_wk.df(d/_t)/(_t*_t) ) : Scalar(0.);
}
