/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
  \file test/common/testUtils.h
  \brief Useful functions for tests
*/

#pragma once

#include "Eigen/Eigen"
#include "Ponca/src/Common/defines.h"
#include "Ponca/src/Common/pointTypes.h"
#include "Ponca/src/Common/pointGeneration.h"
#include PONCA_MULTIARCH_INCLUDE_CU_STD(cmath)

#include <vector>

// Epsilon precision
template<typename T> inline T testEpsilon()
{
    return Eigen::NumTraits<T>::dummy_precision();
}

template<> inline float testEpsilon<float>()
{
    return 1e-2f;
}

template<> inline double testEpsilon<double>()
{
    return 1e-5;
}

template<> inline long double testEpsilon<long double>()
{
    return 1e-5;
}

template<typename DataPoint>
void reverseNormals(std::vector<DataPoint>& _dest, const std::vector<DataPoint>& _src, bool _bRandom = true)
{
    typedef typename DataPoint::VectorType VectorType;

    VectorType vNormal;

    for(unsigned int i = 0; i < _src.size(); ++i)
    {
        vNormal = _src[i].normal();

        if(_bRandom)
        {
            float reverse = Eigen::internal::random<float>(0.f, 1.f);
            if(reverse > 0.5f)
            {
                vNormal = -vNormal;
            }
        }
        else
        {
            vNormal = -vNormal;
        }

        _dest[i] = DataPoint(_src[i].pos(), vNormal);
    }
}

template<typename DataPoint>
typename DataPoint::Scalar getPointKappaMean(typename DataPoint::VectorType _vPoint, typename DataPoint::Scalar _a, typename DataPoint::Scalar _b)
{
    typedef typename DataPoint::Scalar Scalar;

    Scalar x = _vPoint.x();
    Scalar y = _vPoint.y();

    Scalar ax2 = _a * x * _a * x;
    Scalar by2 = _b * y * _b * y;

    Scalar num = (1 + ax2) * _b + (1 + by2) * _a;
    Scalar den = (1 + ax2 + by2);
    den = std::pow(den, Scalar(3./2.));

    Scalar kappa = num / den * Scalar(0.5);

    return kappa;
}

template<typename DataPoint>
typename DataPoint::Scalar getKappaMean(const std::vector<DataPoint>& _vectorPoints, typename DataPoint::VectorType _vCenter,
                                        typename DataPoint::Scalar _a, typename DataPoint::Scalar _b, typename DataPoint::Scalar _analysisScale)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    int size = static_cast<int>(_vectorPoints.size());
    Scalar kappaMean = 0.;
    int nbNei = 0;

    for(int i = 0; i < size; ++i)
    {
        VectorType q = _vectorPoints[i].pos() - _vCenter;

        if(q.norm() <= _analysisScale)
        {
            kappaMean += getPointKappaMean<DataPoint>(_vectorPoints[i].pos(), _a, _b);
            ++nbNei;
        }
    }

    return kappaMean / nbNei;
}

template<typename Fit1, typename Fit2>
void isSamePlane(const Fit1& fit1, const Fit2& fit2) {
    const auto &plane1 = fit1.compactPlane();
    const auto &plane2 = fit2.compactPlane();

    // Test we fit the same plane
    VERIFY(plane1.isApprox(plane2));
}

template<typename Fit1, typename Fit2>
void isSameSphere(const Fit1& fit1, const Fit2& fit2) {
    const auto &sphere1 = fit1.algebraicSphere();
    const auto &sphere2 = fit2.algebraicSphere();

    // Test we fit the same plane
    VERIFY(sphere1 == sphere2);
}

template<typename Fit1, typename Fit2>
void hasSamePlaneDerivatives(const Fit1& fit1, const Fit2& fit2) {
    // Get covariance
    const auto& dpot1 = fit1.covariancePlaneDer().dPotential();
    const auto& dpot2 = fit2.covariancePlaneDer().dPotential();
    const auto& dnor1 = fit1.covariancePlaneDer().dNormal();
    const auto& dnor2 = fit2.covariancePlaneDer().dNormal();

    // Test we compute the same derivatives
    VERIFY(dpot1.isApprox( dpot2 ));
    VERIFY(dnor1.isApprox( dnor2 ));
}
