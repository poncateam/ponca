/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file test/Grenaille/gls_sphere_der.cpp
 \brief Test validity of GLS derivatives for a sphere

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

	typedef typename Fit::VectorArray VectorArray;
    typedef typename Fit::ScalarArray ScalarArray;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(10, 1000);

	Scalar radius = Eigen::internal::random<Scalar>(1,10);
	VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = 10. * std::sqrt(4 * M_PI * radius * radius / nbPoints);

	Scalar epsilon = testEpsilon<Scalar>() * radius;

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center);
    }

	// Test for each point if the Derivatives are equal to 0
    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());

        for(typename vector<DataPoint>::iterator it = vectorPoints.begin();
            it != vectorPoints.end();
            ++it)
        {
            fit.addNeighbor(*it);
        }

        fit.finalize();

		Scalar kappa = fit.kappa();
		ScalarArray dkappa = fit.dkappa();

		VERIFY( (dkappa.array().abs() < epsilon).all() );
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPosistionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
	typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer> FitSmoothOriented;
	typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer> FitConstantOriented;

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

	callSubTests<float, 2>();
	callSubTests<float, 3>();
	callSubTests<double, 3>();
}
