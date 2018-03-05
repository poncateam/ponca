

/*!
    \file test/Grenaille/curvature.cpp
    \brief Test validity curvature estimator
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

    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnPlane<DataPoint>(center,
                                                     direction,
                                                     width,
                                                     _bAddPositionNoise,
                                                     _bAddNormalNoise);
    }

    // Test for each point if TODO(thib)...
#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {
        epsilon = testEpsilon<Scalar>();
        if ( _bAddPositionNoise) // relax a bit the testing threshold
          epsilon = Scalar(0.001*MAX_NOISE);

        Fit fit;
        fit.setWeightFunc(WeightFunc(analysisScale));
        fit.init(vectorPoints[i].pos());
        fit.compute(vectorPoints.cbegin(), vectorPoints.cend());

        if( fit.isStable() ){

            // Check if principal curvature values are null
            VERIFY(fit.k1() < epsilon);
            VERIFY(fit.k2() < epsilon);

            // Check if principal curvature directions lie on the plane
            VERIFY(fit.k1Direction().dot(direction) < epsilon);
            VERIFY(fit.k2Direction().dot(direction) < epsilon);
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

//    typedef Basket<Point, WeightSmoothFunc,   CompactPlane, CovariancePlaneFit, NormalCovarianceCurvature> FitSmoothNormalCovariance;
//    typedef Basket<Point, WeightConstantFunc, CompactPlane, CovariancePlaneFit, NormalCovarianceCurvature> FitConstantNormalCovariance;
    typedef Basket<Point, WeightSmoothFunc,   CompactPlane, CovariancePlaneFit, ProjectedNormalCovarianceCurvature> FitSmoothProjectedNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CompactPlane, CovariancePlaneFit, ProjectedNormalCovarianceCurvature> FitConstantProjectedNormalCovariance;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
//        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance, WeightSmoothFunc>() ));
//        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance, WeightConstantFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothProjectedNormalCovariance, WeightSmoothFunc>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantProjectedNormalCovariance, WeightConstantFunc>() ));
    }
    cout << "Ok..." << endl;
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
