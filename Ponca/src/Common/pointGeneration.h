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
#include "./defines.h"

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

    /// \brief Sample points on a circle in the xy plane
    template<typename VectorType>
    [[nodiscard]] VectorType getPointOnCircle(
        typename VectorType::Scalar _radius,
        VectorType _vCenter
    ) {
        using Scalar = typename VectorType::Scalar;
        // Generate random angle
        double theta = Eigen::internal::random<Scalar>(0,2.*EIGEN_PI);

        // Generate random radius (adjusted for uniform area distribution)
        double r = _radius * std::sqrt(Eigen::internal::random<Scalar>(0,1));

        // Convert to Cartesian coordinates
        VectorType p = _vCenter;
        p.x() += r * std::cos(theta);
        p.y() += r * std::sin(theta);

        return p;
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

        /// Generate z value using the equation z = ax^2 + by^2 + cxy + dx + ey + f
        template<typename Scalar>
        [[nodiscard]] inline Scalar
        getParaboloidZ(Scalar _x,
                       Scalar _y,
                       Scalar _a,
                       Scalar _b,
                       Scalar _c,
                       Scalar _d,
                       Scalar _e,
                       Scalar _f)
        {
            return _a*_x*_x + _b*_y*_y + _c*_x*_y + _d*_x + _e*_y + _f;
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
        template<typename DataPoint>
        inline typename std::enable_if<Ponca::hasNormal<DataPoint>::value, void>::type
        getParaboloidNormal(DataPoint& in,
                    typename DataPoint::Scalar _a,
                    typename DataPoint::Scalar _b,
                    typename DataPoint::Scalar _c,
                    typename DataPoint::Scalar _d,
                    typename DataPoint::Scalar _e,
                    typename DataPoint::Scalar _f)
        {
            static constexpr typename DataPoint::Scalar two {2};
            auto& pos = in.pos();
            in.normal() = typename DataPoint::VectorType(two*_a * pos.x() + _c*pos.y() + _d,
                                                         two*_b * pos.y() + _c*pos.x() + _d,
                                                         -1.).normalized();
        }

        template<typename DataPoint>
        inline typename std::enable_if<! Ponca::hasNormal<DataPoint>::value, void>::type
        getParaboloidNormal(DataPoint& in,
                            typename DataPoint::Scalar _a,
                            typename DataPoint::Scalar _b,
                            typename DataPoint::Scalar _c,
                            typename DataPoint::Scalar _d,
                            typename DataPoint::Scalar _e,
                            typename DataPoint::Scalar _f)
        { }
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

    /// Generate point samples on the primitive z = ax^2 + by^2
    /// Points (x,y) are generated in the interval [-_s, -s]^2
    /// See getParaboloidZ and getParaboloidNormal
    template<typename DataPoint>
    [[nodiscard]] DataPoint
    getPointOnParaboloid(typename DataPoint::Scalar _a,
                         typename DataPoint::Scalar _b,
                         typename DataPoint::Scalar _c,
                         typename DataPoint::Scalar _d,
                         typename DataPoint::Scalar _e,
                         typename DataPoint::Scalar _f,
                         typename DataPoint::Scalar _s,
                         bool _bAddNoise = true)
    {
        typedef typename DataPoint::Scalar Scalar;
        typedef typename DataPoint::VectorType VectorType;

        DataPoint out;

        // generate random position in polar coordinates to get points on the circle
        out.pos() = getPointOnCircle(_s, VectorType({0,0,0}));
        out.pos().z() = internal::getParaboloidZ(out.pos().x(), out.pos().y(), _a, _b, _c, _d, _e, _f);

        // does nothing if the point type does not have a normal field.
        internal::getParaboloidNormal(out, _a, _b, _c, _d, _e, _f);

        if(_bAddNoise) //spherical noise
        {
            out.pos() += VectorType::Random().normalized() * Eigen::internal::random<Scalar>(Scalar(0), Scalar(1. - MIN_NOISE));
        }

        return out;
    }


    template<typename DataPoint, typename Params>
    [[nodiscard]] DataPoint
    getPointOnParaboloid(const Params &_params,
                         typename DataPoint::Scalar _s,
                         bool _bAddNoise = true)
    {
        return getPointOnParaboloid<DataPoint>(_params(0), _params(1), _params(2),
                                               _params(3), _params(4), _params(5), _s, _bAddNoise);
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
