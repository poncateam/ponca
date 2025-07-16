/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/

#pragma once
#include "weightFunc.h"

template<class DataPoint, bool _local>
typename NeighborhoodFrameBase<DataPoint, _local>::VectorType
NeighborhoodFrameBase<DataPoint, _local>::convertToGlobalBasisInternal(const VectorType& _q) const
{
    return _q + m_p;
}

template<class DataPoint, bool _local>
typename NeighborhoodFrameBase<DataPoint, _local>::VectorType
NeighborhoodFrameBase<DataPoint, _local>::convertToLocalBasisInternal(const VectorType& _q) const
{
    return _q - m_p;
}

template<class DataPoint, bool _local>
typename NeighborhoodFrameBase<DataPoint, _local>::VectorType
NeighborhoodFrameBase<DataPoint, _local>::convertToGlobalBasis(const VectorType& _q) const
{
    return _local ? convertToGlobalBasisInternal(_q) : _q;
}

template<class DataPoint, bool _local>
typename NeighborhoodFrameBase<DataPoint, _local>::VectorType
NeighborhoodFrameBase<DataPoint, _local>::convertToLocalBasis(const VectorType& _q) const
{
    return _local ? convertToLocalBasisInternal(_q) : _q;
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::WeightReturnType
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::operator()( const DataPoint& _q) const
{
    Scalar d  = NeighborhoodFrame::convertToLocalBasisInternal(_q.pos()).norm();
    return { (d <= m_t) ? m_wk.f(d/m_t) : Scalar(0.), NeighborhoodFrame::convertToLocalBasis(_q.pos()) };
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::VectorType
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::spacedw( const VectorType& _q,
						                                                  const DataPoint&) const
{
    static_assert(WeightKernel::isDValid, "First order derivatives are required");
    VectorType result = VectorType::Zero();
    VectorType q = NeighborhoodFrame::convertToLocalBasis(_q);
    Scalar d = NeighborhoodFrame::convertToLocalBasisInternal(_q).norm();
    if (d <= m_t && d != Scalar(0.)) result = (q / (d * m_t)) * m_wk.df(d/m_t);
    return result;
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::MatrixType
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::spaced2w( const VectorType& _q,
                                                                           const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    MatrixType result = MatrixType::Zero();
    VectorType q = NeighborhoodFrame::convertToLocalBasis(_q);
    Scalar d = NeighborhoodFrame::convertToLocalBasisInternal(_q).norm();
    if (d <= m_t && d != Scalar(0.))
    {
        Scalar der = m_wk.df(d/m_t);
        result = q*q.transpose()/d*(m_wk.ddf(d/m_t)/m_t - der/d);
        result.diagonal().array() += der;
        result *= Scalar(1.)/(m_t*d);
    }
    return result;
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::Scalar
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::scaledw( const VectorType& _q,
                                                                          const DataPoint&) const
{
    static_assert(WeightKernel::isDValid, "First order derivatives are required");
    Scalar d  = NeighborhoodFrame::convertToLocalBasisInternal(_q).norm();
    return (d <= m_t) ? Scalar( - d*m_wk.df(d/m_t)/(m_t*m_t) ) : Scalar(0.);
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::Scalar
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::scaled2w( const VectorType& _q,
                                                                           const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    Scalar d  = NeighborhoodFrame::convertToLocalBasisInternal(_q).norm();
    return (d <= m_t) ? Scalar(Scalar(2.)*d/(m_t*m_t*m_t)*m_wk.df(d/m_t) +
                               d*d/(m_t*m_t*m_t*m_t)*m_wk.ddf(d/m_t)) :
                        Scalar(0.);
}

template <class DataPoint, class WeightKernel, bool _CenterCoordinates>
typename DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::VectorType
DistWeightFuncBase<DataPoint, WeightKernel, _CenterCoordinates>::scaleSpaced2w( const VectorType& _q,
                                                                                const DataPoint&) const
{
    static_assert(WeightKernel::isDDValid, "Second order derivatives are required");
    VectorType result = VectorType::Zero();
    VectorType q = NeighborhoodFrame::convertToLocalBasis(_q);
    Scalar d = NeighborhoodFrame::convertToLocalBasisInternal(_q).norm();
    if (d <= m_t && d != Scalar(0.))
        result = -q/(m_t*m_t)*(m_wk.df(d/m_t)/d + m_wk.ddf(d/m_t)/m_t);
    return result;
}


