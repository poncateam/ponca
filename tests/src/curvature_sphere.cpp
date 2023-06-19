

/*!
    \file test/Grenaille/curvature_sphere.cpp
    \brief Test validity of curvature estimator for sphere
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
    int nbPoints = Eigen::internal::random<int>(10000, 50000);

    Scalar radius = Eigen::internal::random<Scalar>(1,10);
//    Scalar curvature = Scalar(1.)/radius;

    VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);

    Scalar analysisScale = Eigen::internal::random<Scalar>(0.3, std::sqrt(2.f)) * radius;

    Scalar epsilon = testEpsilon<Scalar>();

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise);
    }

    // Test for each point if principal curvature values are equal to the radius
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); i+=10)
    {
        epsilon = testEpsilon<Scalar>();
        if ( _bAddPositionNoise) // relax a bit the testing threshold
          epsilon = Scalar(0.001*MAX_NOISE);

        epsilon = 0.25; // large threshold

        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints);

        if( fit.isStable() )
        {
            // Check if principal curvature values are equal to the inverse radius
//            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmin()-curvature), Scalar(1.), epsilon) );
//            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmax()-curvature), Scalar(1.), epsilon) );

            // Check if principal curvature directions are tangent to the sphere
            VectorType normal = (vectorPoints[i].pos()-center).normalized();

            if(!Eigen::internal::isMuchSmallerThan(std::abs(fit.kminDirection().dot(normal)), Scalar(1.), epsilon) ||
               !Eigen::internal::isMuchSmallerThan(std::abs(fit.kmaxDirection().dot(normal)), Scalar(1.), epsilon))
            {
                Scalar dot1 = std::abs(fit.kminDirection().dot(normal));
                Scalar dot2 = std::abs(fit.kmaxDirection().dot(normal));
                cout << "dot1 = " << dot1 << endl;
                cout << "dot2 = " << dot2 << endl;
            }

            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kminDirection().dot(normal)), Scalar(1.), epsilon) );
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fit.kmaxDirection().dot(normal)), Scalar(1.), epsilon) );
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

    cout << "Testing with perfect sphere..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothProjectedNormalCovariance, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantProjectedNormalCovariance, WeightConstantFunc>() ));
    }
    cout << "Ok..." << endl;

//    cout << "Testing with noisy sphere..." << endl;
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
