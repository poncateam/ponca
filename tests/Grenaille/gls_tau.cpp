/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/gls_tau.cpp
 \brief Test validity GLS tau param

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

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(10, 1000);
    vector<DataPoint> vectorPoints(nbPoints);
    
	//Random plane parameters
	Scalar centerScale =  Eigen::internal::random<Scalar>(0, 10000);
	VectorType vCenter = VectorType::Random() * centerScale;
	VectorType vPlaneNormal = VectorType::Random().normalized();

    Scalar analysisScale = 100.;
	Scalar epsilon = testEpsilon<Scalar>();

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
		Scalar radius = Eigen::internal::random<Scalar>(-100, 100);
		vectorPoints[i] = getPointOnPlane<DataPoint>(vCenter, vPlaneNormal, radius);
    }

	// Test for each point if the point moved from distance d correspond to tau
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
		// Take a random distance to the plane, not too large to have few points in weightning analysis
		Scalar distanceToPlane = Eigen::internal::random<Scalar>(-25, 25);
		VectorType vEvaluationPoint = vectorPoints[i].pos() + distanceToPlane * vPlaneNormal;

        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
		fit.init(vEvaluationPoint);

        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }

        fit.finalize();

		Scalar fitTau = fit.tau();
		fitTau = abs(fitTau);
		distanceToPlane = abs(distanceToPlane);

		// Test Tau
		VERIFY( abs(distanceToPlane - fitTau) < epsilon);
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPosistionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
	typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam> FitSmoothOriented;
	typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam> FitConstantOriented;
	typedef Basket<Point, WeightSmoothFunc, UnorientedSphereFit, GLSParam> FitSmoothUnoriented;
	typedef Basket<Point, WeightConstantFunc, UnorientedSphereFit, GLSParam> FitConstantUnoriented;

	for(int i = 0; i < g_repeat; ++i)
    {
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
		CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
		CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>() ));
		CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test GLS tau param coherance..." << endl;

	callSubTests<float, 3>();
	callSubTests<double, 3>();
}
