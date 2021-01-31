/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/projection.cpp
 \brief Test validity of the direct projection on an algebraic sphere
 \authors Thibault Lejemble
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeNearestPointQuery.h>

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <chrono>
#include <Eigen/Dense>

//#include <PCP/Math/VectorArray.h>

#include <chrono>
using namespace Eigen;

using namespace std;

using namespace pca;
using namespace Ponca;

using uint = unsigned int;
using Scalar = float;


//Dans le fichier de test de thibault has_duplicate sert a vérifier les données d'entrée pour empecher les duplicatas
// et retourne un erreur si c'est le cas 

template<typename DataPoint, typename WeightFunc>
void testFunction()
{
	// Define related structure
	typedef typename DataPoint::Scalar Scalar;
	typedef typename DataPoint::VectorType VectorType;
	typedef typename DataPoint::QuaternionType QuaternionType;
	typedef Basket<DataPoint, WeightFunc, OrientedSphereFit> Fit;

	//generate samples
	int nbPoints = Eigen::internal::random<int>(100, 1000);
	int nbPointsFit = 50;

	// equal probability of having a plane or a random quadric
	VectorType coeff = 5 * VectorType::Random();
	if (Eigen::internal::random<Scalar>(0., 1.) < Scalar(0.5))
	{
		coeff = VectorType::Zero();
	}
	Scalar width = Eigen::internal::random<Scalar>(1., 10.);
	VectorType center = 1000 * VectorType::Random();

	Scalar zmax = std::abs((coeff[0] + coeff[1]) * width*width);
	Scalar analysisScale = std::sqrt(zmax*zmax + width * width);

	Scalar epsilon = Scalar(10.)*testEpsilon<Scalar>();

	Fit fit;
	fit.setWeightFunc(WeightFunc(analysisScale));
	fit.init(center);

	for (int i = 0; i < nbPointsFit; ++i)
	{
		DataPoint p = getPointOnParaboloid<DataPoint>(VectorType(),     // center (not used)
			coeff,
			QuaternionType(), // (not used)
			width,
			false);           // noise
		p.pos() += center;

		fit.addNeighbor(p);
	}

	fit.finalize();

	if (fit.isStable())
	{
		std::vector<VectorType> samples(nbPoints);
		for (int i = 0; i < nbPoints; ++i)
		{
			VectorType p = center + analysisScale * VectorType::Random();
			samples[i] = p;
			VectorType proj = fit.project(p);

			// check that the projected point is on the surface
			VERIFY(fit.potential(proj) < epsilon);
		}

		auto start1 = std::chrono::system_clock::now();
		for (const auto& p : samples)
			fit.project(p);
		auto end1 = std::chrono::system_clock::now();

		auto start2 = std::chrono::system_clock::now();
		for (const auto& p : samples)
			fit.projectDescent(p);
		auto end2 = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
		std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
		std::cout << "Default: " << elapsed_seconds1.count() << " Descent: " << elapsed_seconds2.count() << "s\n";
		VERIFY(elapsed_seconds1 <= elapsed_seconds2);
	}


}

template<typename Scalar, int Dim>
void callSubTests()
{
	typedef PointPositionNormal<Scalar, Dim> Point;

	typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
	typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

	cout << "Testing with parabola..." << endl;
	for (int i = 0; i < g_repeat; ++i)
	{
		CALL_SUBTEST((testFunction<Point, WeightSmoothFunc>()));
		CALL_SUBTEST((testFunction<Point, WeightConstantFunc>()));
	}
	cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }


	for (int i = 0; i < 10; i++) {

		cout << i % 2 * 5 << endl;
	}


    cout << "Test projection for different baskets..." << endl;

	/*callSubTests<float, 3>();
	callSubTests<double, 3>();
	callSubTests<long double, 3>();*/
	typedef Eigen::Matrix <Scalar, 3, 1> Vector3;
	KdTreeNearestPointQuery<Vector3> q;
}
