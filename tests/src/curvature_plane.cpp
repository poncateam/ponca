

/*!
    \file test/Grenaille/curvature_plane.cpp
    \brief Test validity of curvature estimator for plane
 */

#define FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>

#include "Ponca/src/Fitting/defines.h"
#include "Ponca/src/Fitting/basket.h"
#include "Ponca/src/Fitting/covariancePlaneFit.h"
#include "Ponca/src/Fitting/curvature.h"
#include "Ponca/src/Fitting/mongePatch.h"
#include "Ponca/src/Fitting/weightFunc.h"
#include "Ponca/src/Fitting/weightKernel.h"

#include <Ponca/SpatialPartitioning>

using namespace std;
using namespace Ponca;


template<bool hasFundamentalForms>
struct FundamentalFormTester {
    template<typename Fit>
    static inline void test(const Fit &, typename Fit::Scalar) {}
};

template<>
template<typename Fit>
void FundamentalFormTester<true>::test(const Fit &fit, typename Fit::Scalar epsilon) {
    VERIFY(std::abs(fit.kMeanFromWeingartenMap()) < epsilon);
    VERIFY(std::abs(fit.GaussianCurvatureFromWeingartenMap()) <= epsilon);
}






template<typename DataPoint, typename Fit, typename RefFit>
void testFunction(bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    constexpr bool hasFundamentalForms = hasFirstFundamentalForm<Fit>::value;

    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled plane
    int nbPoints = Eigen::internal::random<int>(1000, 5000);

    // use large width to reduce relative influence of the positional noise
    Scalar width  = Eigen::internal::random<Scalar>(100., 200.);
    Scalar height = width;

    Scalar analysisScale = width; //use a large scale to avoir fitting errors
    Scalar centerScale   = Eigen::internal::random<Scalar>(1,10000);
    VectorType center    = VectorType::Zero(); //Random() * centerScale;

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

    KdTreeDense<DataPoint> tree(vectorPoints);

    // Test for several points if principal curvature values are null
    // The points are generated to not be on the border of the generated point cloud
    // to avoid instabilities issues
    int nbTestPoints = QUICK_TESTS ? 1 :nbPoints/10;

//#pragma omp parallel for
    for(int i = 0; i < nbTestPoints; ++i)
    {
        VectorType queryPos = getPointOnPlane<DataPoint>(center,
                                                         direction,
                                                         width/Scalar(2), // avoid to test on borders
                                                         _bAddPositionNoise,
                                                         _bAddNormalNoise).pos();

        epsilon = testEpsilon<Scalar>();
        if ( _bAddPositionNoise) // relax a bit the testing threshold
          epsilon = Scalar(0.01*MAX_NOISE);

        typename Fit::NeighborFilter nf (queryPos, analysisScale);

        Fit fit;
        fit.setNeighborFilter(nf);
        fit.computeWithIds(tree.range_neighbors(queryPos,analysisScale),vectorPoints);

        RefFit refFit;
        refFit.setNeighborFilter(nf);
        refFit.computeWithIds(tree.range_neighbors(queryPos,analysisScale),vectorPoints);

        // use temporaries to help debugging
        VectorType res    = fit.primitiveGradient(queryPos).normalized();
        Scalar resDot     = res.dot(direction);
        VectorType resRef = refFit.primitiveGradient(queryPos).normalized();
        Scalar resDotRef  = resRef.dot(direction);

        if(Scalar(1.) - std::abs( resDot ) > epsilon){
            VectorType res2 = fit.primitiveGradient(queryPos).normalized();
            std::cout << res2.transpose() << std::endl;
        }

        if( refFit.isStable() ){
            VERIFY(fit.isStable());

            VERIFY(Scalar(1.) - std::abs( resDot ) <= epsilon);
            VERIFY(Scalar(1.) - std::abs( resDotRef ) <= epsilon);

            // use temporaries to help debugging
            Scalar meanCurvature     = fit.kMean();
            Scalar gaussianCurvature = fit.GaussianCurvature();

            // Check if we have a plane
            VERIFY(std::abs(meanCurvature) < epsilon);
            VERIFY(std::abs(gaussianCurvature) <= epsilon);

            FundamentalFormTester<hasFundamentalForms>::test(fit, epsilon);
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
    typedef Basket<Point, WeightSmoothFunc  , MongePatchQuadraticFit> FitCovSmooth;
    typedef Basket<Point, WeightConstantFunc, MongePatchQuadraticFit> FitCovConstant;

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothNormalCovariance, FitSmoothNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantNormalCovariance, FitConstantNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitCovSmooth,FitSmoothNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitCovConstant, FitConstantNormalCovariance>(true) ));
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
