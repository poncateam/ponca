/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file examples/cpp/ponca_basic_cpu.cpp
 * \brief Basic use of Ponca
 * \author: Nicolas Mellado, Gautier Ciaudo
 */

#include <algorithm>
#include <iostream>
#include <vector>

#include <Ponca/Instantiate>
#include <Ponca/Ponca>

using namespace std;
using namespace Ponca;

using MyPoint    = PointPositionNormal<double, 3>;
using Scalar     = MyPoint::Scalar;
using VectorType = MyPoint::VectorType;

// Define related structure
using WeightFunc = DistWeightFunc<MyPoint, SmoothWeightKernel<Scalar>>;
using Fit1       = Basket<MyPoint, WeightFunc, OrientedSphereFit, GLSParam>;
using Fit2       = Basket<MyPoint, WeightFunc, UnorientedSphereFit, GLSParam>;
using Fit3       = BasketDiff<Fit1, FitSpaceDer, OrientedSphereDer, GLSDer, CurvatureEstimatorDer,
                              NormalDerivativeWeingartenEstimator, WeingartenCurvatureEstimatorDer>;
using Fit4       = Basket<MyPoint, WeightFunc, SphereFit, GLSParam>;

template <typename Fit>
void test_fit(Fit& _fit, const KdTree<MyPoint>& tree, const VectorType& _p)
{
    constexpr Scalar tmax = Scalar(100.0);

    // Set a weighting function instance
    _fit.setNeighborFilter({_p, tmax});

    _fit.init();

    // Iterate over samples and _fit the primitive
    for (const int i : tree.rangeNeighbors(_p, tmax))
    {
        _fit.addNeighbor(tree.points()[i]);
    }

    // finalize fitting
    _fit.finalize();

    // Test if the fitting ended without errors
    if (_fit.isStable())
    {
        cout << "Center: [" << _fit.center().transpose() << "] ;  radius: " << _fit.radius() << endl;

        cout << "Pratt normalization" << (_fit.applyPrattNorm() ? " is now done." : " has already been applied.")
             << endl;

        // Play with fitting output
        cout << "Value of the scalar field at the initial point: " << _p.transpose() << " is equal to "
             << _fit.potential(_p) << endl;

        cout << "It's gradient at this place is equal to: " << _fit.primitiveGradient(_p).transpose() << endl;

        cout << "Fitted Sphere: " << endl
             << "\t Tau  : " << _fit.tau() << endl
             << "\t Eta  : " << _fit.eta().transpose() << endl
             << "\t Kappa: " << _fit.kappa() << endl;

        cout << "The initial point " << _p.transpose() << endl
             << "Is projected at   " << _fit.project(_p).transpose() << endl;
    }
}

int main()
{
    // set evaluation point and scale
    VectorType p = VectorType::Random();

    // init input data
    constexpr int n = 10000;
    vector<MyPoint> vecs(n);
    std::generate(vecs.begin(), vecs.end(), getRandomPoint<MyPoint>);

    p = vecs.at(0).pos();

    const KdTreeDense<MyPoint> tree{vecs};

    std::cout << "====================\nOrientedSphereFit:\n";
    Fit1 fit1;
    test_fit(fit1, tree, p);

    std::cout << "\n\n====================\nUnorientedSphereFit:\n";
    Fit2 fit2;
    test_fit(fit2, tree, p);

    std::cout << "\n\n====================\nUnorientedSphereFit:\n";
    Fit3 fit3;
    test_fit(fit3, tree, p);

    if (fit3.isStable())
    {
        cout << "eigen values: " << endl;
        cout << fit3.kmin() << endl;
        cout << fit3.kmax() << endl;
        cout << "eigen vectors: " << endl;
        cout << fit3.kminDirection() << endl << endl;
        cout << fit3.kmaxDirection() << endl;
    }

    std::cout << "\n\n====================\nSphereFit:\n";
    Fit4 fit4;
    test_fit(fit4, tree, p);

    return 0;
}
