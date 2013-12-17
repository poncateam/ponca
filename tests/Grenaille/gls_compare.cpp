/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/gls_compare.cpp
 \brief Test coherance of gls compareTo

 \authors: Gautier Ciaudo
 */


#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
void testFunction(bool bUnoriented = false, bool bAddPositionNoise = false, bool bAddNormalNoise = false)
{
    // Define related structure
	typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

	Scalar radius1 = Eigen::internal::random<Scalar>(1, 5);
	Scalar radius2 = Eigen::internal::random<Scalar>(10, 50);

    Scalar analysisScale = 10.f * std::sqrt( 4.f * M_PI * radius2 * radius2 / nbPoints);

	VectorType center = VectorType::Zero();

	Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> sphere1(nbPoints);
	vector<DataPoint> sphere2(nbPoints);

    for(unsigned int i = 0; i < nbPoints; ++i)
    {
        sphere1[i] = getPointOnSphere<DataPoint>(radius1, center, bAddPositionNoise, bAddNormalNoise, bUnoriented);
		sphere2[i] = getPointOnSphere<DataPoint>(radius2, center, bAddPositionNoise, bAddNormalNoise, bUnoriented);
    }

	// Test for each point if the fitted sphere correspond to the theorical sphere
    for(unsigned int i = 0; i < nbPoints - 1; ++i)
    {
        Fit fit1, fit2, fit3;
        fit1.setWeightFunc(WeightFunc(analysisScale));
		fit2.setWeightFunc(WeightFunc(analysisScale));
		fit3.setWeightFunc(WeightFunc(analysisScale));
        fit1.init(sphere1[i].pos());
		fit2.init(sphere1[i+1].pos());
		fit3.init(sphere2[i].pos());

        for(typename vector<DataPoint>::iterator it = sphere1.begin();
            it != sphere1.end();
            ++it)
        {
            fit1.addNeighbor(*it);
			fit2.addNeighbor(*it);
        }
		for(typename vector<DataPoint>::iterator it = sphere2.begin();
			it != sphere2.end();
			++it)
		{
			fit3.addNeighbor(*it);
		}

		fit1.finalize();
		fit2.finalize();
		fit3.finalize();

		if(fit1.isReady() && fit2.isReady() && fit3.isReady())
		{
			Scalar value1 = fit1.compareTo(fit2);
			Scalar value2 = fit1.compareTo(fit3);

			VERIFY( Eigen::internal::isMuchSmallerThan(value1, 1., epsilon) );
			VERIFY( value2 > epsilon );
		}
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

	cout << "Testing with perfect spheres (oriented / unoriented)..." << endl;
	for(int i = 0; i < g_repeat; ++i)
    {
		//Test with perfect sphere
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
		CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>() ));
		//CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true) ));
		//CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true) ));
	}
	cout << "Ok!" << endl;

	cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
	for(int i = 0; i < g_repeat; ++i)
	{
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(false, true, false) ));
		CALL_SUBTEST(( testFunction<Point, FitConstantOriented, WeightConstantFunc>(false, true, false) ));
		//CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, WeightSmoothFunc>(true, true, true) ));
		//CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, WeightConstantFunc>(true, true, true) ));
    }
	cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test GLS compareTo()..." << endl;

	callSubTests<float, 2>();
	callSubTests<float, 3>();
	callSubTests<double, 3>();
	callSubTests<long double, 2>();
	callSubTests<long double, 3>();
}
