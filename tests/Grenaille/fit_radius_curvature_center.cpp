/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file tets/Grenaille/params_derivatives.cpp
 \brief Test validity of different params

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
	typedef typename DataPoint Point;
	typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(10, 1000);
    
	Scalar radiusScale = Eigen::internal::random<Scalar>(1,10);
	Scalar radius = Eigen::internal::random<Scalar>(0,1) * radiusScale;
    Scalar analysisScale = 10.f * std::sqrt( 4.f * M_PI * radius * radius / nbPoints);
	Scalar centerScale = Eigen::internal::random<Scalar>(1,10);
	VectorType center = VectorType::Random() * centerScale;

	Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();
	// epsilon is relative to the radius size
	Scalar radiusEpsilon = epsilon * radius;


    vector<Point> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<Point>(radius, center);
    }

	// Test for each point if the fitted sphere correspond to the theorical sphere
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());

        for(typename vector<Point>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }

        fit.finalize();

        Scalar fitRadiusKappa = Scalar(abs(1.f / fit.kappa()));
		Scalar fitRadiusAlgebraic = fit.radius();
        typename Point::VectorType fitCenter = fit.center();

		Scalar radiusMax = radius * MAX_NOISE;
		Scalar radiusMin = radius * MIN_NOISE;
		
		// Test procedure
		VERIFY( ((fitCenter - center).array().abs() < (radiusMax - radius) + radiusEpsilon).all() );
		VERIFY( (radiusMin - radiusEpsilon < fitRadiusAlgebraic) && (fitRadiusAlgebraic < radiusMax + radiusEpsilon) );
		// Test reparametrization
        VERIFY( (radiusMin - radiusEpsilon < fitRadiusKappa) && (fitRadiusKappa < radiusMax + radiusEpsilon) );
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

    cout << "Test sphere fitting (radius / center) and GLS curvature for different baskets..." << endl;

	callSubTests<float, 2>();
	callSubTests<float, 3>();
	callSubTests<double, 3>();
}
