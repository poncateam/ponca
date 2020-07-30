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
#include "Ponca/Core"

#include <vector>

#define MIN_NOISE 0.99
#define MAX_NOISE 1.01

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

// Basic point
template<typename _Scalar, int _Dim>
class PointPositionNormal
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim,   1, Eigen::DontAlign>		VectorType;
    typedef Eigen::Matrix<Scalar, Dim+1, 1, Eigen::DontAlign>		HVectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign>	MatrixType;
    typedef Eigen::Quaternion<Scalar, Eigen::DontAlign>			QuaternionType;

    PONCA_MULTIARCH inline PointPositionNormal(  const VectorType &pos = VectorType::Zero(),
                                            const VectorType& normal = VectorType::Zero()
                                        )
        : m_pos(pos), m_normal(normal) {}

    PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    PONCA_MULTIARCH inline const VectorType& normal() const { return m_normal; }

    PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }
    PONCA_MULTIARCH inline VectorType& normal() { return m_normal; }

private:
    VectorType m_pos, m_normal;
};

// Basic point
template<typename _Scalar, int _Dim>
class PointPosition
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim,   1, Eigen::DontAlign>		VectorType;
    typedef Eigen::Matrix<Scalar, Dim+1, 1, Eigen::DontAlign>		HVectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign>	MatrixType;
    typedef Eigen::Quaternion<Scalar, Eigen::DontAlign>			QuaternionType;

    PONCA_MULTIARCH inline PointPosition(  const VectorType &pos = VectorType::Zero() )
        : m_pos(pos) {}

    PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }

    PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }

private:
    VectorType m_pos;
};

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
DataPoint getPointOnSphere(typename DataPoint::Scalar _radius, typename DataPoint::VectorType _vCenter, bool _bAddPositionNoise = true,
                           bool _bAddNormalNoise = true, bool _bReverseNormals = false)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    VectorType vNormal = VectorType::Random().normalized();

    VectorType vPosition = _vCenter + vNormal * _radius; // * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);

    if(_bAddPositionNoise)
    {
        //vPosition = _vCenter + vNormal * _radius * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
        vPosition = vPosition + VectorType::Random().normalized() *
                Eigen::internal::random<Scalar>(Scalar(0.), Scalar(1. - MIN_NOISE));
        vNormal = (vPosition - _vCenter).normalized();
    }

    if(_bAddNormalNoise)
    {
        VectorType vTempPos =  vPosition + VectorType::Random().normalized() *
                Eigen::internal::random<Scalar>(Scalar(0.), Scalar(1. - MIN_NOISE));
        vNormal = (vTempPos - _vCenter).normalized();
    }
    if(_bReverseNormals)
    {
        float reverse = Eigen::internal::random<float>(0.f, 1.f);
        if(reverse > 0.5f)
        {
            vNormal = -vNormal;
        }
    }

    return DataPoint(vPosition, vNormal);
}

/*! \brief Generate points on a plane without normals */
template<typename DataPoint>
DataPoint getPointOnRectangularPlane(
        const typename DataPoint::VectorType& _vPosition,
        const typename DataPoint::VectorType& /*_vNormal*/,
        const typename DataPoint::Scalar& _width,
        const typename DataPoint::Scalar& _height,
        const typename DataPoint::VectorType& _localxAxis,
        const typename DataPoint::VectorType& _localyAxis,
        bool _bAddPositionNoise = true)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    Scalar u = Eigen::internal::random<Scalar>(-_width/Scalar(2),
                                                _width/Scalar(2));
    Scalar v = Eigen::internal::random<Scalar>(-_height/Scalar(2),
                                                _height/Scalar(2));

    VectorType vRandomPosition = _vPosition + u*_localxAxis + v*_localyAxis;

    if(_bAddPositionNoise)
    {
        vRandomPosition = vRandomPosition +
                VectorType::Random().normalized() *
                Eigen::internal::random<Scalar>(0., 1. - MIN_NOISE);
    }

    return DataPoint(vRandomPosition);
}

template<typename DataPoint>
DataPoint getPointOnPlane(typename DataPoint::VectorType _vPosition,
                          typename DataPoint::VectorType _vNormal,
                          typename DataPoint::Scalar _radius,
                          bool _bAddPositionNoise = true,
                          bool _bAddNormalNoise = true,
                          bool _bReverseNormals = false	)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::QuaternionType QuaternionType;

    VectorType vRandom;
    VectorType vRandomDirection = VectorType::Zero();
    VectorType vRandomPoint = VectorType::Zero();
    VectorType vLocalUp = _vNormal;

    do
    {
        vRandom = VectorType::Random().normalized(); // Direction in the unit sphere
        vRandomDirection = vRandom.cross(vLocalUp);
    }
    while(vRandomDirection == VectorType::Zero());

    vRandomDirection = vRandomDirection.normalized();
    vRandomPoint = vRandomDirection * _radius;
    vRandomPoint += _vPosition;

    if(_bAddPositionNoise)
    {
        vRandomPoint = vRandomPoint + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(0., 1. - MIN_NOISE);
    }

    if(_bAddNormalNoise)
    {
        VectorType vLocalLeft = vLocalUp.cross(vRandomDirection);
        VectorType vLocalFront = vLocalLeft.cross(vLocalUp);

        Scalar rotationAngle = Eigen::internal::random<Scalar>(Scalar(-M_PI / 16.), Scalar(M_PI / 16.));
        VectorType vRotationAxis = vLocalLeft;
        QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
        qRotation = qRotation.normalized();
        vLocalUp = qRotation * vLocalUp;

        rotationAngle = Eigen::internal::random<Scalar>(Scalar(-M_PI / 16.), Scalar(M_PI / 16.));
        vRotationAxis = vLocalFront;
        qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
        qRotation = qRotation.normalized();
        vLocalUp = qRotation * vLocalUp;
    }

    if(_bReverseNormals)
    {
        float reverse = Eigen::internal::random<float>(0.f, 1.f);
        if(reverse > 0.5f)
        {
            vLocalUp = -vLocalUp;
        }
    }

    return DataPoint(vRandomPoint, vLocalUp);
}

template<typename Scalar>
inline Scalar
getParaboloidZ(Scalar _x, Scalar _y, Scalar _a, Scalar _b)
{
    return _a*_x*_x + _b*_y*_y;
}

template<typename DataPoint>
DataPoint getPointOnParaboloid(typename DataPoint::VectorType /*_vCenter*/,
                               typename DataPoint::VectorType _vCoef,
                               typename DataPoint::QuaternionType /*_qRotation*/,
                               typename DataPoint::Scalar _analysisScale,
                               bool _bAddNoise = true)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    VectorType vNormal;
    VectorType vPosition;

    Scalar a = _vCoef.x();
    Scalar b = _vCoef.y();
    Scalar x, y, z;

    x = Eigen::internal::random<Scalar>(-_analysisScale, _analysisScale);
    y = Eigen::internal::random<Scalar>(-_analysisScale, _analysisScale);
    z = getParaboloidZ(x, y, a, b);

    vNormal = VectorType((a * x), (b * y), -1.).normalized();

    vPosition.x() = x;
    vPosition.y() = y;
    vPosition.z() = z;

    if(_bAddNoise)
    {
        //spherical noise
        vPosition = vPosition + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(Scalar(0), Scalar(1. - MIN_NOISE));
    }

    //vPosition = _qRotation * vPosition + _vCenter;
    //vNormal = _qRotation * vNormal;

    return DataPoint(vPosition, vNormal);
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

    int size = _vectorPoints.size();
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
