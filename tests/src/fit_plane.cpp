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

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/meanPlaneFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/SpatialPartitioning>

#include <vector>

using namespace std;
using namespace Ponca;

template <bool check>
struct CheckSurfaceVariation {
    template <typename Fit, typename Scalar>
    static inline void run(const Fit& fit, Scalar epsilon){
        VERIFY(fit.surfaceVariation() < epsilon);
    }
};

template <>
template <typename Fit, typename Scalar>
void CheckSurfaceVariation<false>::run(const Fit& /*fit*/, Scalar /*epsilon*/){ }


template<typename DataPoint, typename Fit, bool _cSurfVar>
void testFunction(bool _bUnoriented = false, bool _bAddPositionNoise = false, bool _bAddNormalNoise = false, bool conflictAnnounced = false)
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
      epsilon = Scalar(0.01*MAX_NOISE);
    // Test for each point if the fitted plane correspond to the theoretical plane

    KdTreeDense<DataPoint> tree(vectorPoints);

    // Quick testing is requested for coverage
    int size = QUICK_TESTS ? 1 : int(vectorPoints.size());

#ifdef DEBUG
#pragma omp parallel for
#endif
    for(int i = 0; i < size; ++i)
    {

        Fit fit;
        fit.setNeighborFilter({vectorPoints[i].pos(), analysisScale});
        fit.computeWithIds(tree.rangeNeighbors(vectorPoints[i].pos(),analysisScale),vectorPoints);

        auto ret = fit.getCurrentState();

        switch (ret){
            case STABLE: {
                // Check if the plane orientation is equal to the generation direction
                auto primGrad = fit.primitiveGradient(vectorPoints[i].pos());
                auto gen_dir = std::abs(primGrad.dot(direction));

                VERIFY(Scalar(1.) - gen_dir <= epsilon);

                // Check if the surface variation is small
                CheckSurfaceVariation<_cSurfVar>::run(fit, _bAddPositionNoise ? epsilon*Scalar(10.): epsilon);

                // Check if the query point is on the plane
                if(!_bAddPositionNoise)
                    VERIFY(std::abs(fit.potential(vectorPoints[i].pos())) <= epsilon);
                break;
            }
            case CONFLICT_ERROR_FOUND:
                VERIFY(conflictAnnounced); // raise error only if conflict is detected but not announced
                break;
            default:
                VERIFY(MULTIPASS_PLANE_FITTING_FAILED);
        };
        if(conflictAnnounced) VERIFY(ret == CONFLICT_ERROR_FOUND); //check if we did not detect the conflict
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef Ponca::DistWeightFunc<Point, Ponca::ConstantWeightKernel<Scalar> > WeightConstantFuncLocal;
    typedef Ponca::NoWeightFuncGlobal<Point> NoWeightFuncGlobal;
    typedef Ponca::NoWeightFunc<Point> NoWeightFunc;

    typedef Basket<Point, WeightSmoothFunc, CovariancePlaneFit> CovFitSmooth;
    typedef Basket<Point, WeightConstantFuncLocal, CovariancePlaneFit> CovFitConstant;
    typedef Basket<Point, NoWeightFunc, CovariancePlaneFit> CovFitConstant2;
    typedef Basket<Point, NoWeightFuncGlobal, CovariancePlaneFit> CovFitConstantGlobal;

    typedef Basket<Point, WeightSmoothFunc, MeanPlaneFit> MeanFitSmooth;
    typedef Basket<Point, WeightConstantFuncLocal, MeanPlaneFit> MeanFitConstant;
    typedef Basket<Point, NoWeightFunc, MeanPlaneFit> MeanFitConstant2;
    typedef Basket<Point, NoWeightFuncGlobal, MeanPlaneFit> MeanFitConstantGlobal;

    // test if conflicts are detected
    //! [Conflicting type]
    typedef Basket<Point, NoWeightFuncGlobal, Plane, PlaneFrame,
            MeanNormal, MeanPosition, MeanPlaneFitImpl,
            CovarianceFitBase, CovariancePlaneFitImpl> Hybrid1; //test conflict detection in one direction
    //! [Conflicting type]
    typedef Basket<Point, NoWeightFuncGlobal, Plane, PlaneFrame,
            MeanPosition, CovarianceFitBase, CovariancePlaneFitImpl,
            MeanNormal, MeanPlaneFitImpl> Hybrid2;  //test conflict detection in the second direction

    cout << "Testing with perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect plane
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, true>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, true>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant2, true>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitConstantGlobal, false>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitSmooth, false>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitConstant, false>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitConstant2, false>() ));
        CALL_SUBTEST(( testFunction<Point, MeanFitConstantGlobal, false>() ));
        // Check if fitting conflict is detected
        CALL_SUBTEST(( testFunction<Point, Hybrid1, false>(false, false, false, true) ));
        CALL_SUBTEST(( testFunction<Point, Hybrid2, false>(false, false, false, true) ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position" << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, true>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, true>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant2, true>(false, true, true) ));
        // CALL_SUBTEST(( testFunction<Point, CovFitConstantGlobal, WeightConstantFuncGlobal, true>(false, true, true) ));
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
