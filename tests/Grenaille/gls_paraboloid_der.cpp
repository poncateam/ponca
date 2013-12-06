/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/gls_paraboloid_der.cpp
 \brief Test validity of GLS derivatives for a paraboloid

 \authors: Gautier Ciaudo
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction()
{
    // Define related structure
	typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
	typedef typename DataPoint::QuaternionType QuaternionType;

	typedef typename Fit::VectorArray VectorArray;
    typedef typename Fit::ScalarArray ScalarArray;

    //generate sampled paraboloid
    int nbPoints = Eigen::internal::random<int>(10, 1000);

	Scalar coefScale = Eigen::internal::random<Scalar>(1,10);
    Scalar analysisScale = 10. * std::sqrt(4 * M_PI * coefScale * coefScale / nbPoints);

	VectorType vCenter = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);
	VectorType vCoef = VectorType::Random() * coefScale;
	vCoef.y() = vCoef.x();

	Scalar rotationAngle = Eigen::internal::random<Scalar>(0, 2 * M_PI);
	VectorType vRotationAxis = VectorType::Random().normalized();
	QuaternionType qRotation = QuaternionType(Eigen::AngleAxis<Scalar>(rotationAngle, vRotationAxis));
	qRotation = qRotation.normalized();

	Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
		vectorPoints[i] = getPointOnParaboloid<DataPoint>(vCenter, vCoef, qRotation);
    }

	// Test for each point if the Derivatives are equal to 0
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
		fit.init(VectorType(0., 0., 0.));

        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }

        fit.finalize();

		Scalar tau = fit.tau();
		VectorType eta = fit.eta();
		Scalar kappa = fit.kappa();
		ScalarArray dkappa = fit.dkappa();

		Scalar kappa1 = fit.GLSk1();
		Scalar kappa2 = fit.GLSk2();
		Scalar gaussian = fit.GLSGaussianCurvature();

		Scalar a = vCoef.x();
		Scalar b = vCoef.y();
		Scalar a2 = a * a;
		Scalar b2 = b * b;
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPosistionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
	typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

	typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper> FitSmoothOriented;
	typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper> FitConstantOriented;

	for(int i = 0; i < g_repeat; ++i)
    {
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
		CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test sphere derivatives with GLSParam and OrientedSphereFit..." << endl;

	callSubTests<float, 3>();
	callSubTests<double, 3>();
}
