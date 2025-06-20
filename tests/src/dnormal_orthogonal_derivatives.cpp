/*
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
  \file tests/src/dnormal_orthogonal_derivatives.cpp
  \brief Test validity of the spatial differentiation of normal
*/

#include "../split_test_helper.h"
#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/covariancePlaneFit.h>
#include <Ponca/src/Fitting/mlsSphereFitDer.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace Ponca;

// This class defines the input data format
template<typename _Scalar, int _Dim>
class MyPoint
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1, Eigen::DontAlign>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign> MatrixType;

    PONCA_MULTIARCH inline MyPoint(   const VectorType &pos    = VectorType::Zero(),
                                const VectorType& normal = VectorType::Zero())
        : m_pos(pos), m_normal(normal) {}

    PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    PONCA_MULTIARCH inline const VectorType& normal() const { return m_normal; }

    PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }
    PONCA_MULTIARCH inline VectorType& normal() { return m_normal; }

    static inline MyPoint Random(Scalar radius)
    {
        VectorType scale;
        scale.setRandom();
        VectorType n = VectorType::Random().normalized();
        VectorType p =   scale.asDiagonal() * n * radius   // create an ellipse
                         * Eigen::internal::random<Scalar>(Scalar(0.99),Scalar(1.01)); // add high frequency noise

        n = (scale.asDiagonal().inverse() * n).normalized();

        return MyPoint (p, n);
    };


private:
    VectorType m_pos, m_normal;
};


int kmax =
#ifdef NDEBUG
        g_repeat;
#else
        1;
#endif


/// Helper class used to test fits that are restricted to a specific number of dimensions
template<int Dim>
struct Helper{
    template<typename Scalar> void testRestrictedFits(){}
    template<typename Scalar, typename V, typename M>
    void testCov(Scalar curvature, Scalar epsilon, const V& normal, const M& dNormal){}
};

template<>
template<typename Scalar, typename V, typename M>
void Helper<3>::testCov(Scalar curvature, Scalar epsilon, const V& normal, const M& dN) {
    Eigen::SelfAdjointEigenSolver<M> solver(dN);

    VERIFY(std::abs(solver.eigenvalues()[0]) < epsilon);
    VERIFY(std::abs(solver.eigenvalues()[1] - curvature) < epsilon);
    VERIFY(std::abs(solver.eigenvalues()[2] - curvature) < epsilon);

    VERIFY(std::abs(1 - std::abs(normal.dot(solver.eigenvectors().col(0)))) < epsilon);
    VERIFY(std::abs(normal.dot(solver.eigenvectors().col(1))) < epsilon);
    VERIFY(std::abs(normal.dot(solver.eigenvectors().col(2))) < epsilon);
}


template<typename FitType, typename Functor>
void test_orthoDerivatives(Functor f, bool skipCov = false)
{
    using Point = typename FitType::DataPoint;
    using WeightFunc = typename FitType::WeightFunction;
    using VectorType = typename Point::VectorType;
    using MatrixType = typename Point::MatrixType;
    using Scalar = typename FitType::Scalar;

    //generate samples on a sphere
    int nbPoints = Eigen::internal::random<int>(1000, 3000);

    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    Scalar curvature = Scalar(1.)/radius;
    VectorType center = VectorType::Random() * Eigen::internal::random<Scalar>(1, 10000);
    Scalar analysisScale = Eigen::internal::random<Scalar>(Scalar(0.3), std::sqrt(Scalar(2))) * radius;
    Scalar epsilon = testEpsilon<Scalar>();

    vector<Point> vecs(nbPoints);
    for(unsigned int i = 0; i < vecs.size(); ++i)
        vecs[i] = getPointOnSphere<Point>(radius, center, false, false);

    FitType fit;

    int slice = 10;
#ifdef NDEBUG
#pragma omp parallel for private(fit)
#endif
    for(int k=0; k<int(vecs.size())/slice; ++k)
    {
        fit.setWeightFunc(WeightFunc(vecs[k*slice].pos(), analysisScale));
        fit.init();
        fit.compute(vecs);

        if(fit.isStable())
        {
            auto res = f(fit);
            typename Point::VectorType normal  = res.first;  //fit.primitiveGradient();
            typename Point::MatrixType dN = res.second.template middleCols<Point::Dim>(fit.isScaleDer() ? 1: 0);

            // check that we have unitary normal vector
            VERIFY( normal.norm() - Scalar(1) < epsilon );

            auto proj = (normal.transpose() * dN).eval();
            VERIFY( (proj.array() < epsilon).all() );

            if( ! skipCov)
                Helper<Point::Dim>().testCov(curvature, epsilon, normal, dN);
        }
    }
}


template<>
template<typename Scalar>
void Helper<3>::testRestrictedFits() {
    typedef MyPoint<Scalar,3> Point;
    typedef DistWeightFunc<Point,SmoothWeightKernel<Scalar> > WeightFunc;
    using PlaneFit    = BasketDiff<
            Basket<Point, WeightFunc, CovariancePlaneFit>,
            FitScaleSpaceDer, CovariancePlaneDer>;

    for(int k=0; k<kmax; ++k) {
        CALL_SUBTEST(test_orthoDerivatives<PlaneFit>([](auto&fit){
            return std::make_pair(fit.primitiveGradient(), fit.dNormal());}, true));
    }
};


template<typename Scalar, int Dim>
void _testAdimensionalFits(){
    cout << "Test in dimension " << Dim << std::endl;

    typedef MyPoint<Scalar,Dim> Point;
    typedef DistWeightFunc<Point,SmoothWeightKernel<Scalar> > WeightFunc;
    using SphereFit    = BasketDiff<
            Basket<Point, WeightFunc, OrientedSphereFit>,
            FitScaleSpaceDer, OrientedSphereDer>;
    using MlsSphereFit    = BasketDiff<
            Basket<Point, WeightFunc, OrientedSphereFit>,
            FitScaleSpaceDer, OrientedSphereDer, MlsSphereFitDer>;

    for(int k=0; k<kmax; ++k) {
        CALL_SUBTEST(test_orthoDerivatives<SphereFit>([](auto&fit){
            return std::make_pair(fit.primitiveGradient(), fit.dNormal());}));
        CALL_SUBTEST(test_orthoDerivatives<MlsSphereFit>([](auto&fit){
            return std::make_pair(fit.mlsSphereFitDer().primitiveGradient(),
                                  fit.mlsSphereFitDer().dNormal());}));
    }
}

template<typename Scalar, int Dim>
void testFits(){
    _testAdimensionalFits<Scalar, Dim>();
    Helper<Dim>().template testRestrictedFits<Scalar>();
}


int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test orthogonality between the normal vector and its derivatives..." << endl;

    CALL_SUBTEST_1((testFits<float,  2>()));
    CALL_SUBTEST_2((testFits<double,  2>()));

    CALL_SUBTEST_3((testFits<float,  3>()));
    CALL_SUBTEST_4((testFits<double,  3>()));

    CALL_SUBTEST_5((testFits<float,  4>()));
    CALL_SUBTEST_6((testFits<double,  4>()));

    return EXIT_SUCCESS;
}
