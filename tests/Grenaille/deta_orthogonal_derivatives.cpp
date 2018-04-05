/*
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
  \file tets/Grenaille/deta_orthogonal_derivatives.cpp
  \brief Test validity of the spatial differentiation of eta
*/

#include "../common/testing.h"
#include "Patate/grenaille.h"
#include "Eigen/Eigen"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace Grenaille;

// This class defines the input data format
template<typename _Scalar, int _Dim>
class MyPoint
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1, Eigen::DontAlign>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign> MatrixType;

    MULTIARCH inline MyPoint(   const VectorType &pos    = VectorType::Zero(),
                                const VectorType& normal = VectorType::Zero())
        : m_pos(pos), m_normal(normal) {}

    MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    MULTIARCH inline const VectorType& normal() const { return m_normal; }

    MULTIARCH inline VectorType& pos()    { return m_pos; }
    MULTIARCH inline VectorType& normal() { return m_normal; }

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


template<typename Scalar, int Dim>
void test_orthoEta()
{
    // Define related structure
    typedef MyPoint<Scalar,Dim> Point;
    typedef DistWeightFunc<Point,SmoothWeightKernel<Scalar> > WeightFunc;
    typedef Basket<Point, WeightFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer> Fit;

    Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    // generate sample data
    int n = Eigen::internal::random<int>(10,1000);
    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    Scalar tmax = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius/n);
    vector<Point> vecs (n);

    for(int k=0; k<n; ++k)
    {
        vecs[k] = Point::Random(radius);
    }

    Fit fit;

#pragma omp parallel for private(fit)
    for(int k=0; k<int(vecs.size()); ++k)
    {
        fit.setWeightFunc(WeightFunc(tmax));
        fit.init(vecs[k].pos());
        fit.compute(vecs.cbegin(), vecs.cend());

        if(fit.isStable())
        {
            typename Point::VectorType eta  = fit.eta();
            typename Point::MatrixType deta = fit.deta().template middleCols<Point::Dim>(fit.isScaleDer() ? 1: 0);

            VERIFY( ((eta.transpose() * deta).array() < epsilon).all() );
        }
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
        return EXIT_FAILURE;

    cout << "Test orthogonality between eta and its derivatives..." << endl;

    for(int k=0; k<g_repeat; ++k)
    {
        CALL_SUBTEST(( test_orthoEta<float,  3>() ));
        CALL_SUBTEST(( test_orthoEta<double, 2>() ));
        CALL_SUBTEST(( test_orthoEta<double, 4>() ));
    }

    return EXIT_SUCCESS;
}
