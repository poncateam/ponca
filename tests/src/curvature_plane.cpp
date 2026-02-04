

/*!
    \file test/Grenaille/curvature_plane.cpp
    \brief Test validity of curvature estimator for plane
 */

#define FITTING_FAILED false

#include "../common/testing.h"
#include "../common/testUtils.h"

#include "../split_test_helper.h"

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
    using Scalar = typename Fit::Scalar;

    VERIFY(std::abs(fit.weingartenCurvatureEstimator().kMean()) < epsilon);
    VERIFY(std::abs(fit.weingartenCurvatureEstimator().GaussianCurvature()) <= epsilon);

    VERIFY(std::abs(fit.fundamentalFormWeingartenEstimator().kMean()) < epsilon);
    VERIFY(std::abs(fit.fundamentalFormWeingartenEstimator().GaussianCurvature()) <= epsilon);

    // Check that principal curvature are well computed
    VERIFY((Scalar(.5)*(fit.kmin()+fit.kmax()) - fit.kMean()) < epsilon);
    VERIFY(((fit.kmin()*fit.kmax()) - fit.GaussianCurvature()) < epsilon);
    VERIFY(fit.kmin()< epsilon); // we know that curvatures should be 0
    VERIFY(fit.kmin()< epsilon); // we know that curvatures should be 0

    // Check that principal curvature directions are well-formed
    VERIFY(fit.kminDirection().norm()-Scalar(1) < epsilon);
    VERIFY(fit.kmaxDirection().norm()-Scalar(1) < epsilon);
    VERIFY(fit.kminDirection().dot(fit.kmaxDirection()) < epsilon);
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
        fit.computeWithIds(tree.rangeNeighbors(queryPos,analysisScale),vectorPoints);

        RefFit refFit;
        refFit.setNeighborFilter(nf);
        refFit.computeWithIds(tree.rangeNeighbors(queryPos,analysisScale),vectorPoints);

        // use temporaries to help debugging
        VectorType res    = fit.primitiveGradient(queryPos).normalized();
        Scalar resDot     = res.dot(direction);
        VectorType resRef = refFit.primitiveGradient(queryPos).normalized();
        Scalar resDotRef  = resRef.dot(direction);

        if(Scalar(1.) - std::abs( resDot ) > epsilon){
            VectorType res2 = fit.primitiveGradient(queryPos).normalized();
            // std::cout << res2.transpose() << std::endl;
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

    // Covariance-based fits
    typedef Basket<Point, WeightSmoothFunc  , CovariancePlaneFit> FitSmoothNormalCovariance;
    typedef Basket<Point, WeightConstantFunc, CovariancePlaneFit> FitConstantNormalCovariance;
    // Curvature estimators that runs on top of the covariance fit
    typedef BasketDiff< FitSmoothNormalCovariance, FitSpaceDer, CovariancePlaneDer,
                        CurvatureEstimatorDer, NormalDerivativeWeingartenEstimator> EstimatorSmoothNormalCovariance;
    typedef BasketDiff< FitConstantNormalCovariance, FitSpaceDer, CovariancePlaneDer,
                        CurvatureEstimatorDer, NormalDerivativeWeingartenEstimator> EstimatorConstantNormalCovariance;
    // Curvature estimators based on MongePatch fitting using generalized quadric
    typedef Basket<Point, WeightSmoothFunc  , MongePatchQuadraticFit> FitMongeSmooth;
    typedef Basket<Point, WeightConstantFunc, MongePatchQuadraticFit> FitMongeConstant;
    // Curvature estimators based on MongePatch fitting using restricted quadric
    typedef Basket<Point, WeightSmoothFunc  , MongePatchRestrictedQuadraticFit> FitMongeRestrictedSmooth;
    typedef Basket<Point, WeightConstantFunc, MongePatchRestrictedQuadraticFit> FitMongeRestrictedConstant;

    cout << "Testing with perfect plane..." << flush;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, EstimatorSmoothNormalCovariance, EstimatorSmoothNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, EstimatorConstantNormalCovariance, EstimatorConstantNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitMongeSmooth,FitSmoothNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitMongeConstant, EstimatorConstantNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitMongeRestrictedSmooth,FitSmoothNormalCovariance>() ));
        CALL_SUBTEST(( testFunction<Point, FitMongeRestrictedConstant, EstimatorConstantNormalCovariance>() ));
    }
    cout << "Ok..." << endl;

    cout << "Testing with noisy plane..." << flush;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, EstimatorSmoothNormalCovariance, EstimatorSmoothNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, EstimatorConstantNormalCovariance, EstimatorConstantNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitMongeSmooth,EstimatorSmoothNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitMongeConstant, EstimatorConstantNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitMongeRestrictedSmooth,EstimatorSmoothNormalCovariance>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitMongeRestrictedConstant, EstimatorConstantNormalCovariance>(true) ));
    }
    cout << "Ok..." << endl;
}


int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test curvature estimation with float" << endl;
    CALL_SUBTEST_1((callSubTests<float, 3>()));
    cout << "Test curvature estimation with double" << endl;
    CALL_SUBTEST_2((callSubTests<double, 3>()));
    cout << "Test curvature estimation with long double" << endl;
    CALL_SUBTEST_3((callSubTests<long double, 3>()));
}

