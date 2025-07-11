/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file test/src/fit_cov.cpp
 \brief Test validity of the covariance fitting procedures wrt to standard algorithm
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/defines.h>
#include <Ponca/src/Fitting/mean.h>
#include <Ponca/src/Fitting/covarianceFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <vector>

using namespace std;
using namespace Ponca;


/// Class that perform the covariance fit using standard two-passes procedure
template < class DataPoint, class _WFunctor, typename T>
class CovarianceFitTwoPassesBase : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_MEAN_POSITION,
        PROVIDES_POSITION_COVARIANCE
    };

public:
    using MatrixType = typename DataPoint::MatrixType; /*!< \brief Alias to matrix type*/
    /*! \brief Solver used to analyse the covariance matrix*/
    using Solver = Eigen::SelfAdjointEigenSolver<MatrixType>;

protected:
    // computation data
    MatrixType m_cov {MatrixType::Zero()};     /*!< \brief Covariance matrix */
    Solver m_solver;  /*!<\brief Solver used to analyse the covariance matrix */
    bool m_barycenterReady {false};
    VectorType m_barycenter;
    Scalar sumW;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(CovarianceFitTwoPassesBase,covarianceFitTwoPasses)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    /*! \brief Reading access to the Solver used to analyse the covariance matrix */
    PONCA_MULTIARCH inline const Solver& solver() const { return m_solver; }
};

template < class DataPoint, class _WFunctor, typename T>
void
CovarianceFitTwoPassesBase<DataPoint, _WFunctor, T>::init()
{
    Base::init();
    m_cov.setZero();
    m_barycenterReady = false;
    m_barycenter.setZero();
    sumW = Scalar(0);
}

template < class DataPoint, class _WFunctor, typename T>
bool
CovarianceFitTwoPassesBase<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                             const VectorType &localQ,
                                                             const DataPoint &attributes)
{
    if( ! m_barycenterReady )  /// first pass
    {
        return Base::addLocalNeighbor(w, localQ, attributes);
    }
    else {                     /// second pass
        VectorType q = localQ - m_barycenter; // saved from previous run
        m_cov  += w * q * q.transpose();
        sumW   += w;
        return true;
    }
    return false;
}


template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
CovarianceFitTwoPassesBase<DataPoint, _WFunctor, T>::finalize ()
{
    if( ! m_barycenterReady ){   /// end of the first pass
        auto ret = Base::finalize();
        if(ret == STABLE) {
            m_barycenterReady = true;
            m_barycenter = Base::barycenterLocal();
            return NEED_OTHER_PASS;
        }
        // handle specific configurations
        if( ret != STABLE)
            return Base::m_eCurrentState;
        }
    else {                       /// end of the second pass
        // Normalize covariance matrix
        m_cov /= sumW ;  /// \warning There is a bug in the pass system that prevents to call Base::getWeightSum();

        m_solver.compute(m_cov);
        Base::m_eCurrentState = ( m_solver.info() == Eigen::Success ? STABLE : UNDEFINED );
    }

    return Base::m_eCurrentState;
}


template<typename DataPoint, typename Fit, typename FitRef, typename WeightFunc, bool _cSurfVar> //, typename Fit, typename WeightFunction>
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

#pragma omp parallel for
    for(int i = 0; i < int(vectorPoints.size()); ++i)
    {

        Fit fit;
        fit.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        auto fitState = fit.compute(vectorPoints);

        FitRef ref;
        ref.setWeightFunc(WeightFunc(vectorPoints[i].pos(), analysisScale));
        ref.init();
        auto refState = ref.compute(vectorPoints);

        VERIFY(fitState == refState);

        auto checkVectors = [](const VectorType& x, const VectorType& y){
            Scalar values =  std::min(
                    (x.array() - y.array()).abs().matrix().squaredNorm(),
                    (x.array() + y.array()).abs().matrix().squaredNorm()); // deal with sign ambiguity
            VERIFY( Eigen::internal::isApproxOrLessThan(values, Scalar(1e-5)) );
        };

        // check we get the same decomposition
        checkVectors(fit.covarianceFit().solver().eigenvalues(),
                     ref.covarianceFitTwoPasses().solver().eigenvalues());
        for (int d = 0; d != 3; ++d)
            checkVectors(fit.covarianceFit().solver().eigenvectors().col(d),
                         ref.covarianceFitTwoPasses().solver().eigenvectors().col(d));
        }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPositionNormal<Scalar, Dim> Point;

    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightSmoothFunc;
    typedef DistWeightFunc<Point, ConstantWeightKernel<Scalar> > WeightConstantFunc;

    typedef Basket<Point, WeightSmoothFunc, MeanPosition, CovarianceFitBase> CovFitSmooth;
    typedef Basket<Point, WeightConstantFunc, MeanPosition, CovarianceFitBase> CovFitConstant;

    typedef Basket<Point, WeightSmoothFunc, PrimitiveBase, MeanPosition, CovarianceFitTwoPassesBase> RefFitSmooth;
    typedef Basket<Point, WeightConstantFunc, PrimitiveBase, MeanPosition, CovarianceFitTwoPassesBase> RefFitConstant;


    cout << "Testing with data sampling a perfect plane..." << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        //Test with perfect plane
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, RefFitSmooth, WeightSmoothFunc, true>() ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, RefFitConstant, WeightConstantFunc, true>() ));
    }
    cout << "Ok!" << endl;

    cout << "Testing with noise on position" << endl;
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, CovFitSmooth, RefFitSmooth, WeightSmoothFunc, true>(false, true, true) ));
        CALL_SUBTEST(( testFunction<Point, CovFitConstant, RefFitConstant, WeightConstantFunc, true>(false, true, true) ));
    }
    cout << "Ok!" << endl;
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test covariance matrix construction and decomposition for different baskets..." << endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
}
