/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

/*!
  \file pointGeneration.h
  \brief Point generation methods
*/
#include "Eigen/Eigen"

#define MIN_NOISE 0.99
#define MAX_NOISE 1.01

namespace Ponca {

    /*! \brief Generate points on a sphere */
    template<typename DataPoint>
    [[nodiscard]] DataPoint getPointOnSphere(
        const typename DataPoint::Scalar _radius,
        const typename DataPoint::VectorType _vCenter,
        const bool _bAddPositionNoise = true,
        const bool _bAddNormalNoise = true,
        const bool _bReverseNormals = false
    ) {
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;

        VectorType vNormal   = VectorType::Random().normalized();
        VectorType vPosition = _vCenter + vNormal * _radius;

        if(_bAddPositionNoise)
        {
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
            const float reverse = Eigen::internal::random<float>(0.f, 1.f);
            if(reverse > 0.5f)
                vNormal = -vNormal;
        }

        return DataPoint(vPosition, vNormal);
    }

    /*! \brief Generate points on a plane without normals */
    template<typename DataPoint>
    [[nodiscard]] DataPoint getPointOnRectangularPlane(
        const typename DataPoint::VectorType& _vPosition,
        const typename DataPoint::VectorType& /*_vNormal*/,
        const typename DataPoint::Scalar& _width,
        const typename DataPoint::Scalar& _height,
        const typename DataPoint::VectorType& _localxAxis,
        const typename DataPoint::VectorType& _localyAxis,
        const bool _bAddPositionNoise = true
    ) {
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;

        const Scalar u = Eigen::internal::random<Scalar>(
            -_width /Scalar(2), _width /Scalar(2) );
        const Scalar v = Eigen::internal::random<Scalar>(
            -_height/Scalar(2), _height/Scalar(2) );

        VectorType vRandomPosition = _vPosition + u*_localxAxis + v*_localyAxis;

        if(_bAddPositionNoise)
        {
            vRandomPosition = vRandomPosition +
                    VectorType::Random().normalized() *
                    Eigen::internal::random<Scalar>(0., 1. - MIN_NOISE);
        }

        return DataPoint(vRandomPosition);
    }

    /*! \brief Generate points on a plane */
    template<typename DataPoint>
    [[nodiscard]] DataPoint getPointOnPlane(
        const typename DataPoint::VectorType _vPosition,
        const typename DataPoint::VectorType _vNormal,
        const typename DataPoint::Scalar _radius,
        const bool _bAddPositionNoise = true,
        const bool _bAddNormalNoise = true,
        const bool _bReverseNormals = false
    ) {
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;
        typedef Eigen::Quaternion<Scalar> QuaternionType;

        VectorType vRandom;
        VectorType vRandomDirection = VectorType::Zero();
        VectorType vRandomPoint = VectorType::Zero();
        VectorType vLocalUp     = _vNormal;

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
            vRandomPoint = vRandomPoint + VectorType::Random().normalized() *
                    Eigen::internal::random<Scalar>(Scalar(0), Scalar(1. - MIN_NOISE));
        }

        if(_bAddNormalNoise)
        {
            VectorType vLocalLeft = vLocalUp.cross(vRandomDirection);
            VectorType vLocalFront = vLocalLeft.cross(vLocalUp);

            Scalar rotationAngle     = Eigen::internal::random<Scalar>(Scalar(-M_PI / 16.), Scalar(M_PI / 16.));
            VectorType vRotationAxis = vLocalLeft;
            QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
            qRotation = qRotation.normalized();
            vLocalUp  = qRotation * vLocalUp;

            rotationAngle = Eigen::internal::random<Scalar>(Scalar(-M_PI / 16.), Scalar(M_PI / 16.));
            vRotationAxis = vLocalFront;
            qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
            qRotation = qRotation.normalized();
            vLocalUp = qRotation * vLocalUp;
        }

        if(_bReverseNormals)
        {
            const float reverse = Eigen::internal::random<float>(0.f, 1.f);
            if(reverse > 0.5f)
                vLocalUp = -vLocalUp;

        }

        return DataPoint(vRandomPoint, vLocalUp);
    }

    namespace internal {
        /*! \brief Generate z value using the equation z = ax^2 + by^2 */
        template<typename Scalar>
        [[nodiscard]] inline Scalar getParaboloidZ(
            const Scalar _x,
            const Scalar _y,
            const Scalar _a,
            const Scalar _b
        ) {
            return _a*_x*_x + _b*_y*_y;
        }

        /*! \brief Generate z value using the equation z = ax^2 + by^2 */
        template<typename VectorType>
        [[nodiscard]] inline VectorType getParaboloidNormal(
            const VectorType& in,
            const typename VectorType::Scalar _a,
            const typename VectorType::Scalar _b
        ) {
            return VectorType((_a * in.x()), (_b * in.y()), -1.).normalized();;
        }
    }

    /*! \brief Generate point samples on the primitive z = ax^2 + by^2
     *
     * Points (x,y) are generated in the interval [-_s, -s]^2
     */
    template<typename DataPoint>
    [[nodiscard]] DataPoint getPointOnParaboloid(
        const typename DataPoint::Scalar _a,
        const typename DataPoint::Scalar _b,
        const typename DataPoint::Scalar _s,
        const bool _bAddNoise = true
    ) {
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;

        VectorType vNormal;
        VectorType vPosition;

        const Scalar x = Eigen::internal::random<Scalar>(-_s, _s),
                     y = Eigen::internal::random<Scalar>(-_s, _s);

        vPosition = { x, y, internal::getParaboloidZ(x, y, _a, _b)};
        vNormal   = internal::getParaboloidNormal(vPosition, _a, _b);

        if(_bAddNoise) // spherical noise
            vPosition += VectorType::Random().normalized() * Eigen::internal::random<Scalar>(Scalar(0), Scalar(1. - MIN_NOISE));

        return DataPoint(vPosition, vNormal);
    }

    /*! \brief Generate a random points */
    template<typename DataPoint>
    [[nodiscard]] DataPoint getRandomPoint() {
        typedef typename DataPoint::VectorType VectorType;
        typedef typename DataPoint::Scalar Scalar;
        const VectorType n = VectorType::Random().normalized();
        const VectorType p = n * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
        return {p, (n + VectorType::Random()*Scalar(0.1)).normalized()};
    }
}
