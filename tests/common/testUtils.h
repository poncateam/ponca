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
	return Eigen::NumTraits<Scalar>::dummy_precision();
}

template<> inline float testEpsilon<float>()
{
	return 1e-3f;
}

template<> inline double testEpsilon<double>()
{
	return 1e-3f;
}

// Basic point
template<typename _Scalar, int _Dim>
class PointPosistionNormal
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

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
DataPoint getPointOnSphere(typename DataPoint::Scalar radius, typename DataPoint::VectorType vCenter, bool bAddNoise = true)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    VectorType vNormal = VectorType::Random().normalized();

    VectorType vPosition;
    if(bAddNoise)
    {
        vPosition = vCenter + vNormal * radius * Eigen::internal::random<Scalar>(MIN_NOISE, MAX_NOISE);
    }
	else
	{
		vPosition = vCenter + vNormal * radius;
	}

    //vNormal = vPosition.normalized();

    return DataPoint(vPosition, vNormal);
}

template<typename DataPoint>
DataPoint getPointOnPlane(typename DataPoint::VectorType vPosition, typename DataPoint::VectorType vNormal, typename DataPoint::Scalar radius)
{
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

	VectorType vRandom;
	VectorType vRandomPoint = VectorType::Zero();

	do
	{
		vRandom = VectorType::Random().normalized(); // Direction in the unit sphere
		vRandomPoint = vRandom.cross(vNormal);
	}
	while(vRandomPoint == VectorType::Zero());

	vRandomPoint = vRandomPoint.normalized();
	vRandomPoint *= radius;
	vRandomPoint += vPosition;

    return DataPoint(vRandomPoint, vNormal);
}
#endif // _TEST_UTILS_H_
