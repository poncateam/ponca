

/*!
    \file test/Grenaille/curvature_plane.cpp
    \brief Test validity of curvature estimator for plane
 */

#define FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

using namespace std;
using namespace Grenaille;


template<typename DataPoint, typename Fit, typename WeightFunc> //, typename Fit, typename WeightFunction>
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

    // Test for each point if principal curvature values are null
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        epsilon = testEpsilon<Scalar>();
        if ( _bAddPositionNoise) // relax a bit the testing threshold
          epsilon = Scalar(0.001*MAX_NOISE);

        Fit fit;
        fit.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        fit.init();
        fit.compute(vectorPoints);

        if( fit.isStable() ){

            // Check if principal curvature values are null
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmin()), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmax()), Scalar(1.), epsilon) );

            // Check if principal curvature directions lie on the plane
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kminDirection().dot(direction)), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmaxDirection().dot(direction)), Scalar(1.), epsilon) );
        }
        else {
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

    typedef Basket<Point, WeightSmoothFunc,   CompactPlane, CovariancePlaneFit, NormalCovarianceCurvature> FitSmoothNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CompactPlane, CovariancePlaneFit, NormalCovarianceCurvature> FitConstantNormalCovariance;
    typedef Basket<Point, WeightSmoothFunc,   CompactPlane, CovariancePlaneFit, ProjectedNormalCovarianceCurvature> FitSmoothProjectedNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CompactPlane, CovariancePlaneFit, ProjectedNormalCovarianceCurvature> FitConstantProjectedNormalCovariance;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothProjectedNormalCovariance, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantProjectedNormalCovariance, WeightConstantFunc>() ));
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
