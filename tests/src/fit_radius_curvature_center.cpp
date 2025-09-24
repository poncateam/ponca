/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/Grenaille/fit_radius_curvature_center.cpp
 \brief Test validity of algebraic sphere procedure and GLS kappa
 */


#include "../split_test_helper.h"
#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/gls.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/unorientedSphereFit.h>
#include <Ponca/src/Fitting/curvature.h>
#include <Ponca/src/Fitting/curvatureEstimation.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;

template <bool isSpaceDer, int Dim>
struct subTestSpatial{
    template<typename Scalar, typename Fit>
    static
    inline void eval(const Fit& /*_fit*/,
                     Scalar /*_radiusMin*/,
                     Scalar /*_radiusMax*/,
                     Scalar /*_radiusEpsilon*/,
                     Scalar /*_fitRadiusKappa*/,
                     Scalar /*_fitRadiusAlgebraic*/) { }
};

template <>
template<typename Scalar, typename Fit>
void
subTestSpatial<true, 3>::eval(const Fit& _fit,
                           Scalar _radiusMin,
                           Scalar _radiusMax,
                           Scalar _radiusEpsilon,
                           Scalar /*_fitRadiusKappa*/,
                           Scalar /*_fitRadiusAlgebraic*/){

    Scalar kmin      = _fit.kmin();
    Scalar kmax      = _fit.kmax();
    Scalar kmean   = (kmin + kmax) / Scalar (2.);
    Scalar radius  = Scalar(1.) / kmean;

    //Test value
    VERIFY( (radius > _radiusMin - _radiusEpsilon) &&
            (radius < _radiusMax + _radiusEpsilon) );
}



template<typename DataPoint, typename Fit, bool isSpaceDer>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false)
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;

    //generate sampled sphere
    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(2., 10.);

    Scalar analysisScale = Scalar(10.) * std::sqrt( Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale = Eigen::internal::random<Scalar>(1,10000);
    VectorType center = VectorType::Random() * centerScale;

    Scalar epsilon = testEpsilon<Scalar>();
    // epsilon is relative to the radius size
    Scalar radiusEpsilon = epsilon * radius;


    vector<DataPoint> vectorPoints(nbPoints);

    for(unsigned int i = 0; i < vectorPoints.size(); ++i)
    {
        vectorPoints[i] = getPointOnSphere<DataPoint>(radius, center, _bAddPositionNoise, _bAddNormalNoise, _bUnoriented);
    }

    // Quick testing is requested for coverage
    int size = QUICK_TESTS ? 1 : int(vectorPoints.size());

    // Test for each point if the fitted sphere correspond to the theoretical sphere
#ifdef NDEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < size; ++i)
    {
        Fit fit;
        fit.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit.compute(vectorPoints);

        if(fit.isStable())
        {
            Scalar fitRadiusKappa = Scalar(std::abs(Scalar(1.) / fit.kappa()));
            Scalar fitRadiusAlgebraic = fit.radius();
            VectorType fitCenter (fit.center());

            Scalar radiusMax = radius * Scalar(MAX_NOISE);
            Scalar radiusMin = radius * Scalar(MIN_NOISE);

            // Test procedure
            VERIFY( (fitCenter - center).norm() < (radiusMax - radius) + radiusEpsilon );
            VERIFY( (fitRadiusAlgebraic > radiusMin - radiusEpsilon) && (fitRadiusAlgebraic < radiusMax + radiusEpsilon) );
            // Test re-parameterization
            VERIFY( (fitRadiusKappa > radiusMin - radiusEpsilon) && (fitRadiusKappa < radiusMax + radiusEpsilon) );
            //Test coherence
            VERIFY( Eigen::internal::isMuchSmallerThan(std::abs(fitRadiusAlgebraic - fitRadiusKappa), Scalar(1.), epsilon) );

            //Test using spatial derivatives if defined
            subTestSpatial<isSpaceDer, DataPoint::Dim>::eval(fit, radiusMin, radiusMax, radiusEpsilon, fitRadiusKappa, fitRadiusAlgebraic);

            //Test on eta
            if(!_bAddPositionNoise && !_bAddNormalNoise)
            {
                //sometimes eta can be reversed
                VectorType fitEta (fit.eta().normalized().array().abs());
                VectorType expectedEta (vectorPoints[i].normal().array().abs());

                VERIFY( Eigen::internal::isMuchSmallerThan((fitEta - expectedEta).norm(), Scalar(1.), epsilon)  );
            }
        }
    }
}

#define DECLARE_DEFAULT_TYPES                                                                      \
    typedef PointPositionNormal<Scalar, Dim> Point;                                                \
                                                                                                   \
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;                   \
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;               \
                                                                                                   \
    typedef Basket<Point, WeightSmoothFunc, OrientedSphereFit, GLSParam> FitSmoothOriented;        \
    typedef Basket<Point, WeightConstantFunc, OrientedSphereFit, GLSParam> FitConstantOriented;    \
    typedef Basket<Point, WeightSmoothFunc, UnorientedSphereFit, GLSParam> FitSmoothUnoriented;    \
    typedef Basket<Point, WeightConstantFunc, UnorientedSphereFit, GLSParam> FitConstantUnoriented;

template<typename Scalar, int Dim>
void callSubTests()
{
    DECLARE_DEFAULT_TYPES

    cout << "Testing with perfect sphere (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect sphere
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, false>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, false>() ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, false>(true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, false>(true) ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented)..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOriented, false>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOriented, false>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitSmoothUnoriented, false>(true, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantUnoriented, false>(true, true, true) ));
    }
    cout << "Ok!" << endl;
}

template<typename Scalar, int Dim>
void callDerivativeSubTests()
{
    DECLARE_DEFAULT_TYPES

    using FitSmoothOrientedSpatial   = BasketDiff<FitSmoothOriented, FitSpaceDer, OrientedSphereDer, CurvatureEstimatorBaseDiff, NormalDerivativesCurvatureEstimator>;
    using FitConstantOrientedSpatial = BasketDiff<FitConstantOriented, FitSpaceDer, OrientedSphereDer, CurvatureEstimatorBaseDiff, NormalDerivativesCurvatureEstimator>;

    cout << "Testing with perfect sphere (oriented / unoriented) with spatial derivatives..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect sphere
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpatial, true>() ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpatial, true>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position and normals (oriented / unoriented) with spatial derivatives..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, FitSmoothOrientedSpatial, true>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, FitConstantOrientedSpatial, true>(false, true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test sphere fitting (radius / center) and GLS curvature for different baskets..." << endl;

    CALL_SUBTEST_1(( callSubTests<float,       2>() ));
    CALL_SUBTEST_2(( callSubTests<float,       3>() ));
    CALL_SUBTEST_3(( callSubTests<double,      3>() ));
    CALL_SUBTEST_4(( callSubTests<long double, 2>() ));
    CALL_SUBTEST_5(( callSubTests<long double, 3>() ));

    CALL_SUBTEST_6(( callDerivativeSubTests<float,       3>() ));
    CALL_SUBTEST_7(( callDerivativeSubTests<double,      3>() ));
    CALL_SUBTEST_8(( callDerivativeSubTests<long double, 3>() ));
}
