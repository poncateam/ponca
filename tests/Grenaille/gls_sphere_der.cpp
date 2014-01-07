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
void testFunction(bool bAddPositionNoise = false, bool bAddNormalNoise = false)
{
    // Define related structure
	typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

	typedef typename Fit::VectorArray VectorArray;
    typedef typename Fit::ScalarArray ScalarArray;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

	Scalar radius = Eigen::internal::random<Scalar>(1,10);
	VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = 10. * std::sqrt(4 * M_PI * radius * radius / nbPoints);

	Scalar epsilon = testEpsilon<Scalar>();
	Scalar radisuEpsilon = epsilon * radius;

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, bAddPositionNoise, bAddNormalNoise);
    }

	// Test for each point if the Derivatives of kappa are equal to 0
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

		if(fit.isStable())
		{
			ScalarArray dkappa = fit.dkappa();

			for(int i = 0; i < dkappa.size(); ++i)
			{
				VERIFY( Eigen::internal::isMuchSmallerThan(dkappa[i], 1., epsilon) );
			}
		}
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPosistionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;

    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer> FitSmoothOriented;

	cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
	for(int i = 0; i < g_repeat; ++i)
    {
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>() ));
    }
	cout << "Ok!" << endl;

	/*cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
	for(int i = 0; i < g_repeat; ++i)
	{
		CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, WeightSmoothFunc>(true, true) ));
	}
	cout << "Ok!" << endl;*/
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test sphere derivatives with GLSParam and OrientedSphereFit..." << endl;

	callSubTests<double, 2>();
	callSubTests<float, 3>();
	callSubTests<long double, 3>();
	callSubTests<double, 3>();
}
