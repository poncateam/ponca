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
#include <Ponca/src/Common/pointTypes.h>
#include <Ponca/src/Common/pointGeneration.h>

#include "Eigen/Eigen"

#include <vector>

using namespace std;
using namespace Ponca;

using MyPoint = PointPositionNormal<double, 3>;
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
  std::generate(vecs.begin(), vecs.end(), getRandomPoint<MyPoint>);

  std::cout << "====================\nCovariancePlaneFit:\n";
  CovPlaneFit fit;
  test_fit(fit, vecs, vecs.at(0).pos());
}
