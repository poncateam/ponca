/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::VectorType
DistWeightFunc<DataPoint, WeightKernel>::convertToLocalBasis(const VectorType& _q) const
{
    return _q - m_p;
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::WeightReturnType
DistWeightFunc<DataPoint, WeightKernel>::w( const VectorType& _q, 
					                        const DataPoint&) const
{
    VectorType q = convertToLocalBasis(_q);
    Scalar d  = q.norm();
    return { (d <= m_t) ? m_wk.f(d/m_t) : Scalar(0.), q };
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::VectorType
DistWeightFunc<DataPoint, WeightKernel>::spacedw(   const VectorType& _q, 
						                            const DataPoint&) const
{
    static_assert(WeightKernel::isDValid, "First order derivatives are required");
    VectorType result = VectorType::Zero();
    VectorType q = convertToLocalBasis(_q);
    Scalar d = q.norm();
    if (d <= m_t && d != Scalar(0.)) result = (q / (d * m_t)) * m_wk.df(d/m_t);
    return result;
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::MatrixType
DistWeightFunc<DataPoint, WeightKernel>::spaced2w(   const VectorType& _q,
                                                     const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    MatrixType result = MatrixType::Zero();
    VectorType q = convertToLocalBasis(_q);
    Scalar d = q.norm();
    if (d <= m_t && d != Scalar(0.))
    {
        Scalar der = m_wk.df(d/m_t);
        result = q*q.transpose()/d*(m_wk.ddf(d/m_t)/m_t - der/d);
        result.diagonal().array() += der;
        result *= Scalar(1.)/(m_t*d);
    }
    return result;
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::scaledw(   const VectorType& _q, 
						                            const DataPoint&) const
{
    static_assert(WeightKernel::isDValid, "First order derivatives are required");
    Scalar d  = convertToLocalBasis(_q).norm();
    return (d <= m_t) ? Scalar( - d*m_wk.df(d/m_t)/(m_t*m_t) ) : Scalar(0.);
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::Scalar
DistWeightFunc<DataPoint, WeightKernel>::scaled2w(   const VectorType& _q,
                                                     const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    Scalar d  = convertToLocalBasis(_q).norm();
    return (d <= m_t) ? Scalar(Scalar(2.)*d/(m_t*m_t*m_t)*m_wk.df(d/m_t) +
                               d*d/(m_t*m_t*m_t*m_t)*m_wk.ddf(d/m_t)) :
                        Scalar(0.);
}

template <class DataPoint, class WeightKernel>
typename DistWeightFunc<DataPoint, WeightKernel>::VectorType
DistWeightFunc<DataPoint, WeightKernel>::scaleSpaced2w(   const VectorType& _q,
                                                          const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    VectorType result = VectorType::Zero();
    VectorType q = convertToLocalBasis(_q);
    Scalar d = q.norm();
    if (d <= m_t && d != Scalar(0.))
        result = -q/(m_t*m_t)*(m_wk.df(d/m_t)/d + m_wk.ddf(d/m_t)/m_t);
    return result;
}


