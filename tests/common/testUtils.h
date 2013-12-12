/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _TEST_UTILS_H_
#define _TEST_UTILS_H_

#include "Eigen/Eigen"
#include "Patate/grenaille.h"

#define MIN_NOISE 0.99f
#define MAX_NOISE 1.01f

// Epsilon precision
template<typename T> inline T testEpsilon()
{
	return Eigen::NumTraits<T>::dummy_precision();
}

template<> inline float testEpsilon<float>()
{
	return 1e-3f;
}

template<> inline double testEpsilon<double>()
{
	return 1e-6f;
}

template<> inline long double testEpsilon<long double>()
{
	return 1e-6f;
}

// Basic point
template<typename _Scalar, int _Dim>
class PointPosistionNormal
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1, Eigen::DontAlign>		VectorType;
	typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign>	MatrixType;
	typedef Eigen::Quaternion<Scalar, Eigen::DontAlign>			QuaternionType;

    MULTIARCH inline PointPosistionNormal(   const VectorType &pos = VectorType::Zero(), 
											 const VectorType& normal = VectorType::Zero()
										 )
                                : _pos(pos), _normal(normal) {}

    MULTIARCH inline const VectorType& pos()    const { return _pos; }  
    MULTIARCH inline const VectorType& normal() const { return _normal; }

    MULTIARCH inline VectorType& pos()    { return _pos; }  
    MULTIARCH inline VectorType& normal() { return _normal; }

private:
    VectorType _pos, _normal;
};

template<typename DataPoint>
DataPoint getPointOnSphere(typename DataPoint::Scalar radius, typename DataPoint::VectorType vCenter, bool bAddPositionNoise = true,
						   bool bAddNormalNoise = true, bool bReverseNormals = false)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    VectorType vNormal = VectorType::Random().normalized();

    VectorType vPosition = vCenter + vNormal * radius; // * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);

    if(bAddPositionNoise)
    {
		VectorType vPosition = vCenter + vNormal * radius * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
		//vPosition = vPosition + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
		vNormal = (vPosition - vCenter).normalized();
    }

	if(bAddNormalNoise)
	{
		VectorType vTempPos = vPosition + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
		vNormal = (vTempPos - vCenter).normalized();
	}

	if(bReverseNormals)
	{
		float reverse = Eigen::internal::random<float>(0.f, 1.f);
		if(reverse > 0.5f)
		{
			vNormal = -vNormal;
		}	
	}

    return DataPoint(vPosition, vNormal);
}

template<typename DataPoint>
DataPoint getPointOnPlane(typename DataPoint::VectorType vPosition, typename DataPoint::VectorType vNormal, typename DataPoint::Scalar radius,
						  bool bAddPositionNoise = true, bool bAddNormalNoise = true, bool bReverseNormals = false	)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
	typedef typename DataPoint::QuaternionType QuaternionType;

	VectorType vRandom;
	VectorType vRandomDirection = VectorType::Zero();
	VectorType vRandomPoint = VectorType::Zero();
	VectorType vLocalUp = vNormal;

	do
	{
		vRandom = VectorType::Random().normalized(); // Direction in the unit sphere
		vRandomDirection = vRandom.cross(vLocalUp);
	}
	while(vRandomDirection == VectorType::Zero());

	vRandomDirection = vRandomPoint.normalized();
	vRandomPoint = vRandomDirection * radius;
	vRandomPoint += vPosition;

	if(bAddPositionNoise)
	{
		vRandomPoint = vRandomPoint + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
	}

	if(bAddNormalNoise)
	{
		VectorType vLocalLeft = vLocalUp.cross(vRandomDirection);
		VectorType vLocalFront = vLocalLeft.cross(vLocalUp);

		Scalar rotationAngle = Eigen::internal::random<Scalar>(-M_PI / 16., M_PI / 16.);
		VectorType vRotationAxis = vLocalLeft;
		QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
		qRotation = qRotation.normalized();
		vLocalUp = qRotation * vLocalUp;

		rotationAngle = Eigen::internal::random<Scalar>(-M_PI / 16., M_PI / 16.);
		vRotationAxis = vLocalFront;
		qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
		qRotation = qRotation.normalized();
		vLocalUp = qRotation * vLocalUp;
	}

	if(bReverseNormals)
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
Scalar getParaboloidZ(Scalar x, Scalar y, Scalar a, Scalar b)
{
	Scalar x2 = x * x;
	Scalar y2 = y * y;

	Scalar z = (a * x2 + b * y2); //* Scalar(0.5);

	return z;
}

template<typename DataPoint>
DataPoint getPointOnParaboloid(typename DataPoint::VectorType vCenter, typename DataPoint::VectorType vCoef,
							   typename DataPoint::QuaternionType qRotation, bool bAddNoise = true)
{
	typedef typename DataPoint::Scalar Scalar;
	typedef typename DataPoint::VectorType VectorType;

	VectorType vNormal;
	VectorType vPosition;

	Scalar a = vCoef.x();
	Scalar b = vCoef.y();
	Scalar x, y, z;

	//do
	//{
		x = Eigen::internal::random<Scalar>(-0.0001 * a, 0.0001 * a);
		y = Eigen::internal::random<Scalar>(-0.0001 * b, 0.0001 * b);
		z = getParaboloidZ(x, y, a, b);
	//}
	//while(z > Scalar(10.));

	vNormal = VectorType((2. * a * x), (2. * b * y), 1.).normalized();

	vPosition.x() = x;
	vPosition.y() = y;
	vPosition.z() = z;

	if(bAddNoise)
	{
		//spherical noise
		vPosition = vPosition + VectorType::Random().normalized() * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
	}

	//vPosition = qRotation * vPosition + vCenter;
	//vNormal = qRotation * vNormal;

	return DataPoint(vPosition, vNormal);
}
#endif // _TEST_UTILS_H_
