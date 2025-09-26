

/*!
    \file test/Grenaille/curvature_plane.cpp
    \brief Test validity of curvature estimator for plane
 */

#define FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

#include "Ponca/src/Fitting/basket.h"
#include "Ponca/src/Fitting/covariancePlaneFit.h"
#include "Ponca/src/Fitting/curvature.h"
#include "Ponca/src/Fitting/curvatureEstimation.h"
#include "Ponca/src/Fitting/mongePatch.h"
#include "Ponca/src/Fitting/weightFunc.h"
#include "Ponca/src/Fitting/weightKernel.h"

using namespace std;
using namespace Ponca;


template<typename DataPoint, typename Fit>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(1000, 5000);

    Scalar width  = Eigen::internal::random<Scalar>(1., 10.);
    Scalar height = width;

    Scalar analysisScale = Scalar(15.) * std::sqrt( width * height / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1,10000);
    VectorType center    = VectorType::Random() * centerScale;

    VectorType direction = VectorType::Random().normalized();

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise);
    }

    // Quick testing is requested for coverage
    int size = QUICK_TESTS ? 1 : int(vectorPoints.size());

    // Test for each point if principal curvature values are null
#pragma omp parallel for
    for(int i = 0; i < size; ++i)
    {
        epsilon = testEpsilon<Scalar>();
        const auto& queryPos = vectorPoints[i].pos();
        if ( _bAddPositionNoise) // relax a bit the testing threshold
          epsilon = Scalar(0.001*MAX_NOISE);

        Fit fit;
        fit.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit.compute(vectorPoints);

        if( fit.isStable() ){
            VERIFY(Scalar(1.) - std::abs(fit.primitiveGradient(queryPos).dot(direction)) <= epsilon);
            // Check if we have a plane
            VERIFY(std::abs(fit.kMean()) < epsilon);
            VERIFY(fit.GaussianCurvature() <= epsilon);
        } else {
            VERIFY(FITTING_FAILED);
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc  , CovariancePlaneFit, CurvatureEstimatorBase> FitSmoothNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit, CurvatureEstimatorBase> FitConstantNormalCovariance;
    typedef Basket<Point, WeightSmoothFunc  , CovariancePlaneFit, CurvatureEstimatorBase> FitSmoothProjectedNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit, CurvatureEstimatorBase> FitConstantProjectedNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit, CurvatureEstimatorBase> FitConstantProjectedNormalCovariance;
    typedef Basket<Point, WeightSmoothFunc  , CovariancePlaneFit, MongePatch> FitCovSmooth;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit, MongePatch> FitCovConstant;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothProjectedNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantProjectedNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitCovSmooth>() ));
        CALL_SUBTEST(( testFunction<Point, FitCovConstant>() ));
    }
    cout << "Ok..." << endl;

//    cout << "Testing with noisy plane..." << endl;
//    for(int i = 0; i < g_repeat; ++i)
//    {
//        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance, WeightSmoothFunc>(true, true) ));
//        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance, WeightConstantFunc>(true, true) ));
//        CALL_SUBTEST(( testFunction<Point, FitSmoothProjectedNormalCovariance, WeightSmoothFunc>(true, true) ));
//        CALL_SUBTEST(( testFunction<Point, FitConstantProjectedNormalCovariance, WeightConstantFunc>(true, true) ));
//    }
//    cout << "Ok..." << endl;
}


int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test curvature estimation..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
