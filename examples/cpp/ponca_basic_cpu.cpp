/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
\file examples/Grenaille/basic_cpu.h
\brief Basic use of Grenaille

\author: Nicolas Mellado, Gautier Ciaudo
*/
#include <cmath>
#include <algorithm>
#include <iostream>

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/gls.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/unorientedSphereFit.h>
#include <Ponca/src/Fitting/sphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>
#include <Ponca/src/Fitting/curvatureEstimation.h>
#include <Ponca/src/Fitting/curvature.h>

#include <Ponca/SpatialPartitioning>

#include "Eigen/Eigen"

#include <vector>

using namespace std;
using namespace Ponca;


// This class defines the input data format
class MyPoint
{
public:
  enum {Dim = 3};
  typedef double Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
  typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

  PONCA_MULTIARCH inline MyPoint(const VectorType& _pos    = VectorType::Zero(),
                                 const VectorType& _normal = VectorType::Zero())
    : m_pos(_pos), m_normal(_normal) {}

  PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const VectorType& normal() const { return m_normal; }

  PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }
  PONCA_MULTIARCH inline VectorType& normal() { return m_normal; }

  static inline MyPoint Random()
  {
    VectorType n = VectorType::Random().normalized();
    VectorType p = n * Eigen::internal::random<Scalar>(0.9,1.1);
    return MyPoint (p, (n + VectorType::Random()*0.1).normalized());
  }

private:
  VectorType m_pos, m_normal;
};

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

// Define related structure
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
using Fit1 = Basket<MyPoint,WeightFunc,OrientedSphereFit,   GLSParam>;
using Fit2 = Basket<MyPoint,WeightFunc,UnorientedSphereFit, GLSParam>;
using Fit3 = BasketDiff< Fit1, FitSpaceDer, OrientedSphereDer, GLSDer, CurvatureEstimatorBase, NormalDerivativesCurvatureEstimator>;
using Fit4 = Basket<MyPoint,WeightFunc,SphereFit, GLSParam>;

template<typename Fit>
void test_fit(Fit& _fit, const KdTree<MyPoint>& tree, const VectorType& _p)
{
  Scalar tmax = 100.0;

  // Set a weighting function instance
  _fit.setWeightFunc(WeightFunc(tmax));

  // Set the evaluation position
  _fit.init(_p);

  // Iterate over samples and _fit the primitive
  for(int i : tree.range_neighbors(_p, tmax) )
  {
      _fit.addNeighbor( tree.points()[i] );
  }

  //finalize fitting
  _fit.finalize();

  //Test if the fitting ended without errors
  if(_fit.isStable())
  {
    cout << "Center: [" << _fit.center().transpose() << "] ;  radius: " << _fit.radius() << endl;

    cout << "Pratt normalization"
      << (_fit.applyPrattNorm() ? " is now done." : " has already been applied.") << endl;

    // Play with fitting output
    cout << "Value of the scalar field at the initial point: "
      << _p.transpose()
      << " is equal to " << _fit.potential(_p)
      << endl;

    cout << "It's gradient at this place is equal to: "
      << _fit.primitiveGradient(_p).transpose()
      << endl;

    cout << "Fitted Sphere: " << endl
      << "\t Tau  : "      << _fit.tau()             << endl
      << "\t Eta  : "      << _fit.eta().transpose() << endl
      << "\t Kappa: "      << _fit.kappa()           << endl;

    cout << "The initial point " << _p.transpose()              << endl
      << "Is projected at   " << _fit.project(_p).transpose() << endl;
  }

}

int main()
{
  // set evaluation point and scale
  VectorType p = VectorType::Random();

  // init input data
  int n = 10000;
  vector<MyPoint> vecs (n);
  std::generate(vecs.begin(), vecs.end(), []() {return MyPoint::Random(); });

  p = vecs.at(0).pos();

  KdTree<MyPoint> tree {vecs};

  std::cout << "====================\nOrientedSphereFit:\n";
  Fit1 fit1;
  test_fit(fit1, tree, p);

  std::cout << "\n\n====================\nUnorientedSphereFit:\n";
  Fit2 fit2;
  test_fit(fit2, tree, p);

  std::cout << "\n\n====================\nUnorientedSphereFit:\n";
  Fit3 fit3;
  test_fit(fit3, tree, p);

  if(fit3.isStable())
  {
    cout << "eigen values: "<< endl;
    cout << fit3.kmin() << endl;
    cout << fit3.kmax() << endl;
    cout << "eigen vectors: "<< endl;
    cout << fit3.kminDirection() << endl << endl;
    cout << fit3.kmaxDirection() << endl;
  }

  std::cout << "\n\n====================\nSphereFit:\n";
  Fit4 fit4;
  test_fit(fit4, tree, p);

  return 0;
}
