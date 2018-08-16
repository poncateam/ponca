/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/fit_plane.cpp
 \brief Test validity of plane fitting procedure(s)
 */

#define MULTIPASS_PLANE_FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;

template <bool check>
struct CheckSurfaceVariation {
    template <typename Fit, typename Scalar>
    static inline void run(const Fit& fit, Scalar epsilon){
        VERIFY(fit.surfaceVariation() < epsilon);
    }
};

template <>
template <typename Fit, typename Scalar>
void
CheckSurfaceVariation<false>::run(const Fit& /*fit*/, Scalar /*epsilon*/){ }


template<typename DataPoint, typename Fit, typename WeightFunc, bool _cSurfVar> //, typename Fit, typename WeightFunction>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar height = width;

    Scalar analysisScale = Scalar(15.) * std::sqrt( width * height / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1,10000);
    VectorType center    = VectorType::Random() * centerScale;

    VectorType direction = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();
    // epsilon is relative to the radius size
    //Scalar radiusEpsilon = epsilon * radius;


//    // Create a local basis on the plane by crossing with global XYZ basis
//    // This test can be useful to check if we can retrieve the local basis
//    static const VectorType yAxis(0., 1., 0.);
//    static const VectorType zAxis(0., 0., 1.);
//    VectorType localxAxis =
//            Scalar(1.) - direction.dot(yAxis) > Scalar(0.1) // angle sufficient to cross
//            ? direction.cross(yAxis)
//            : direction.cross(zAxis);
//    VectorType localyAxis = direction.cross(localxAxis);


    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise,
                                                     _bUnoriented);
    }

    epsilon = testEpsilon<Scalar>();
    if ( _bAddPositionNoise) // relax a bit the testing threshold
      epsilon = Scalar(0.001*MAX_NOISE);
    // Test for each point if the fitted plane correspond to the theoretical plane
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {

        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints.cbegin(), vectorPoints.cend());

        if( fit.isStable() ){

            // Check if the plane orientation is equal to the generation direction
            VERIFY(Scalar(1.) - std::abs(fit.primitiveGradient(vectorPoints[i].pos()).dot(direction)) <= epsilon);

            // Check if the surface variation is small
            CheckSurfaceVariation<_cSurfVar>::run(fit, _bAddPositionNoise ? epsilon*Scalar(10.): epsilon);

            // Check if the query point is on the plane
            if(!_bAddPositionNoise)
              VERIFY(fit.potential(vectorPoints[i].pos()) <= epsilon);

        }
        else {
            VERIFY(MULTIPASS_PLANE_FITTING_FAILED);
        }

//        if(fit.isStable())
//        {
//            Scalar fitRadiusKappa = Scalar(std::abs(Scalar(1.) / fit.kappa()));
//            Scalar fitRadiusAlgebraic = fit.radius();
//            VectorType fitCenter = fit.center();

//            Scalar radiusMax = radius * MAX_NOISE;
//            Scalar radiusMin = radius * MIN_NOISE;

//            // Test procedure
//            VERIFY( (fitCenter - center).norm() < (radiusMax - radius) + radiusEpsilon );
//            VERIFY( (fitRadiusAlgebraic > radiusMin - radiusEpsilon) && (fitRadiusAlgebraic < radiusMax + radiusEpsilon) );
//            // Test reparametrization
//            VERIFY( (fitRadiusKappa > radiusMin - radiusEpsilon) && (fitRadiusKappa < radiusMax + radiusEpsilon) );
//            //Test coherance
//            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fitRadiusAlgebraic - fitRadiusKappa), 1., epsilon) );

//            //Test on eta
//            if(!_bAddPositionNoise && !_bAddNormalNoise)
//            {
//                //sometimes eta can be reversed
//                VectorType fitEta = fit.eta().normalized().array().abs();
//                VectorType theoricEta = vectorPoints[i].normal().array().abs();

//                VERIFY( Eigen::internal::isMuchSmallerThan((fitEta - theoricEta).norm(), 1., epsilon)  );
//            }
//        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, CompactPlane, CovariancePlaneFit> CovFitSmooth;
    typedef Basket<Point, WeightConstantFunc, CompactPlane, CovariancePlaneFit> CovFitConstant;

    typedef Basket<Point, WeightSmoothFunc, CompactPlane, MeanPlaneFit> MeanFitSmooth;
    typedef Basket<Point, WeightConstantFunc, CompactPlane, MeanPlaneFit> MeanFitConstant;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect plane
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, WeightSmoothFunc, true>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, WeightConstantFunc, true>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitSmooth, WeightSmoothFunc, false>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitConstant, WeightConstantFunc, false>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position" << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, WeightSmoothFunc, true>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, WeightConstantFunc, true>(false, true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test plane fitting for different baskets..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
