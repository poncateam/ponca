/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
\file examples/Grenaille/basic_cpu.h
\brief Basic use of Grenaille with plane fit

\author: Nicolas Mellado, Gautier Ciaudo
*/
#include <cmath>
#include <algorithm>
#include <iostream>

#include <Ponca/Fitting>

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

  PONCA_MULTIARCH inline MyPoint(   const VectorType& _pos    = VectorType::Zero(),
                                    const VectorType& _normal = VectorType::Zero()
      )
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
  };

private:
  VectorType m_pos, m_normal;
};

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

// Define related structure
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint, WeightFunc, CovariancePlaneFit> CovPlaneFit;

template<typename Fit>
void test_fit(Fit& _fit, vector<MyPoint>& _vecs, const VectorType& _p)
{
  Scalar tmax = 100.0;

  // Set a weighting function instance
  _fit.setNeighborFilter({_p, tmax});

  // Fit plane (method compute handles multipass fitting
  _fit.compute( _vecs.begin(), _vecs.end() );

  if(_fit.isStable())
  {
        cout << "Value of the scalar field at the initial point: "
            << _p.transpose()
            << " is equal to " << _fit.potential(_p)
            << endl;

        cout << "It's gradient at this place is equal to: "
            << _fit.primitiveGradient(_p).transpose()
            << endl;

        cout << "The initial point " << _p.transpose()              << endl
            << "Is projected at   " << _fit.project(_p).transpose() << endl;

        cout << "Value of the surface variation: "
            << _fit.surfaceVariation()
            << endl;
  }
}

int main()
{
  // init random point cloud
  int n = 10000;
  vector<MyPoint> vecs (n);
  std::generate(vecs.begin(), vecs.end(), []() {return MyPoint::Random(); });

  std::cout << "====================\nCovariancePlaneFit:\n";
  CovPlaneFit fit;
  test_fit(fit, vecs, vecs.at(0).pos());
}
